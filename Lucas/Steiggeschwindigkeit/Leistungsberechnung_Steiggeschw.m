%% Leistungsberechnung f�r Flugsysteme

%% Initialisierungen

% Batterie

U_Bat_nom = N_Bat_cell * U_Bat_cell;        % nominale Batteriespannung
U_Bat_min = N_Bat_cell * U_Bat_cell_min;    % minimale Batteriespannung
% C_Bat = E_Dichte * m_Bat / U_Bat_nom;       % Kapazitaet der Batterie in As
C_Bat = N_Bat_cell_p*C_Bat_cell*3600;
% Delta_C_Bat = 0;                            % Initialisierung Batteriekapazit�t, die nach jedem delta_h gebraucht wird



% Propeller

%Entnahme des Durchmessers und des Pitches aus dem Propellernahmen
D = str2double(prop_name(1:strfind(prop_name,'x')-1));      % Durchmesser extrahieren
P_75 = prop_name(strfind(prop_name,'x')+1:end);             % Pitch extrahieren
while isnan(str2double(P_75)) == 1
        P_75(end) = [];
end
P_75 = str2double(P_75);                    % Pitch festlegen
R = D * 0.0254 / 2;                         % Propellerradius in Meter
F = pi * R^2;                               % Fl�che eines Propellers in Quadratmeter
Theta_75 = atan( 4*P_75 / (3*pi * D) );     % geometrischer Anstellwinkel des Propellers bei 75% des Radius
[RPM_map, V_map, T_map, P_map, TAU_map] = Propeller_map(DATA_APC,prop_name);    % Aufbau des Kennfeldes





% Umgebungsparameter

T_11 = T_0 - 0.0065 *(11000-H_0);                   % T in 11000m H�he                     
rho_11 = rho_0 * (1 - 0.0065*(11000/T_0))^4.256;    % Dichte in 11000m H�he
p_11 = p_0 * (1 - 0.0065*(11000/T_0))^5.256;        % Druck in 11000m H�he



% Intialisierung der Matrizen f�r jeden H�henabschnitt

lengthh = floor(abs(H_max - H_0) / Delta_H + 1); 
lengthvkg = floor(abs(V_Kg_min - V_Kg_max) / V_Kg_Delta + 1);
H = zeros(lengthh,1);
U_mot = zeros(lengthh,1);
I_mot = zeros(lengthh,1);
Thrust = zeros(lengthh,1);
Omega = zeros(lengthh,1);
tau = zeros(lengthh,1);
C_Rate = zeros(lengthh,1);
alpha = zeros(lengthh,1);
C_Rest_V = zeros(lengthh,1);
Theta = zeros(lengthh,1);
rho = zeros(lengthh,1);
I_Bat = zeros(lengthh,1);
P_Bat = zeros(lengthh,1);
PWM = zeros(lengthh,1);
M_tip = zeros(lengthh,1);
eta_prop = zeros(lengthh,1);
eta_ges = zeros(lengthh,1);
V_Kg_opt = zeros(lengthh,1);
Delta_C_Bat = zeros(lengthh,1);
Delta_C_Bat_inter = zeros(lengthvkg,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Programmanfang

q = 1;                                      % Z�hler intialisiern
for h_variabel = H_0:Delta_H:H_max
    
    %% Umgebungsparameter als Funktion der H�he
    % Berechnung der Flugh�he f�r iterative Schritte und der mittleren
    % Dichte zwischen den Diskretisierungspunkten
    
    H_unten = h_variabel;
    H_oben = h_variabel + Delta_H;
    H_mitte = (H_oben + H_unten)/2;
    
    % Berechnung der Dichte an den Intevallgrenzen nach Normatmosph�renbedingungen
    
    if H_oben <= 11000
        
        rho_unten = rho_0 *(1-0.0065*(H_unten/T_0))^4.256;
        rho_oben = rho_0 *(1-0.0065*(H_oben/T_0))^4.256;
        
        p_unten = p_0 *(1-0.0065*(H_unten/T_0))^5.256;
        p_oben = p_0 *(1-0.0065*(H_oben/T_0))^5.256;
        
    elseif H_unten <= 11000 && H_oben >11000
        
        rho_unten = rho_0 *(1 - 0.0065*(H_unten/T_0))^4.256;
        rho_oben = rho_11 * exp(-g/(287*T_11)*(H_oben-11000));
        
        p_unten = p_0 *(1-0.0065*(H_unten/T_0))^5.256;
        p_oben = p_11 * exp(-g/(287.1*T_11)*(H_oben-11000));
        
    else
        
        rho_unten = rho_11 * exp(-g/(287*T_11)*(H_unten-11000));
        rho_oben = rho_11 * exp(-g/(287*T_11)*(H_oben-11000));
        
        p_unten = p_11 * exp(-g/(287.1*T_11)*(H_unten-11000));
        p_oben = p_11 * exp(-g/(287.1*T_11)*(H_oben-11000));
        
    end
    
    if H_mitte <= 11000                                             % mittlere Temperatur im Intervall
        T = T_0 - 0.0065 * H_mitte;
    else
        T = T_11;
    end
    
    rho(q) = rho_unten + (rho_oben - rho_unten)/2;                  % Berechnung der mittleren Dichte im Intervall
    p = (p_oben + p_unten)/2;                                       % mittlerer Druck im Intervall
    a = sqrt(kappa*p/rho(q));                                       % Schallgeschwindigkeit in m/s
    
    
    % Dichte an der oberen (_2) und unteren (_1) Intervallgrenze
    if q == 1
        rho_1 = rho_0;
        rho_2 = rho(q);
    else
        rho_1 = rho(q-1);
        rho_2 = rho(q);
    end
    
    
    T_map = T_map * rho(q)/rho_1;                                   % Anpassung des Schubkennfeldes an die sich �ndernde Dichte
    P_map = P_map * rho(q)/rho_1;                                   % Anpassung des Leistungskennfeldes an die sich �ndernde Dichte
    TAU_map = TAU_map * rho(q)/rho_1;                               % Anpassung des Drehmomentkennfeldes an die sich �ndernde Dichte
    
    
    
    
    %% Initialisierungen f�r die Bahngeschwindigkeit

    V_Kg = zeros(lengthvkg,1);
    P_Untersuchung = zeros(lengthvkg,1);
    Thrust_inter = zeros(lengthvkg,1);
    Theta_inter = zeros(lengthvkg,1);
    alpha_inter = zeros(lengthvkg,1);
    Omega_inter = zeros(lengthvkg,1);
    I_mot_inter = zeros(lengthvkg,1);
    U_mot_inter = zeros(lengthvkg,1);
    C_Rate_inter = zeros(lengthvkg,1);
    C_Rest_V_inter = zeros(lengthvkg,1);
    I_Bat_inter = zeros(lengthvkg,1);
    tau_inter = zeros(lengthvkg,1);
    M_tip_inter = zeros(lengthvkg,1);
    PWM_inter = zeros(lengthvkg,1);
    eta_ges_inter = zeros(lengthvkg,1);
    
    %%  Hier beginnt Untersuchung der Bahngeschwindigkeit  
    z = 1;
    
    for V_Kg_variabel = V_Kg_min:V_Kg_Delta:V_Kg_max
        
        % Bahngeschwindigkeit festlegen
        V_Kg(z) = V_Kg_variabel;
        
        % Berechne hier alle weiteren Gr��en des Flugsystems
        t_Flug = Delta_H / V_Kg(z);                                                % Flugzeit
        
        if Abfrage_Flugsystem == 1                                              % handelt es sich um Multicopter oder Fl�chenflugzeug
            
            % MULTICOPTER
            m = m_copter + m_Bat + m_Mot * n_Prop + m_nutz;                     % Gesamtmasse des Quadrocopters
            
            % Aerodynamik berechnen
            
            [Thrust_inter(z),Theta_inter(z),V_A,alpha_inter(z)] = MulticopterAerodynamik(u_Wg,V_Kg(z),gamma,c_W_copter_seitlich,c_W_copter_oben,c_A_copter_max,rho(q),A_copter,m,g);
            
        else
            
            % FL�CHENFLUGZEUG
            m = m_flugzeug + m_Bat + m_Mot * n_Prop + m_nutz;                   % Gesamtmasse des Fl�chenflugzeugs
            
            % Aerodynamik
            
            [Thrust_inter(z),V_A] = FlaechenflugzeugAerodynamik(m,g,epsilon,V_Kg(z));
            
            
        end
        
        
        Thrust_inter(z) = Thrust_inter(z) / n_Prop;                                         % Schub auf n Propeller verteilen
        
        
        
        if Thrust_inter(z) > max(max(T_map))                                          % wenn Schub zu gross (Ergebnis verwerfen)
            Omega_inter(z) = NaN;
            I_mot_inter(z) = NaN;
            C_Rate_inter(z) = NaN;
            C_Rest_V_inter(z) = NaN;
        else
            
            
            % Drehzahl und Drehmoment bestimmen
            
            [Omega_inter(z),tau_inter(z)] = Propeller(V_A, alpha_inter(z), Thrust_inter(z), RPM_map, V_map, T_map, TAU_map);
            
            
            % Wie gro� ist die Blattspitzengeschwindigkeit?
            
            M_tip_inter(z) = (Omega_inter(z) * R)/a;                                        % Blattspitzengeschwindigkeit in Ma
            
            
            % Motorzustand berechnen
            
            [U_mot_inter(z),I_mot_inter(z)] = Motor(tau_inter(z),K_V,I_0,R_i,Omega_inter(z));
            
            
            % Zustand der Motorregler berechnen
            
            [PWM_inter(z),eta_PWM] = ESC(U_mot_inter(z),U_Bat_nom);
            
            
            % Batteriezustand berechnen
            
            [I_Bat_inter(z),C_Rate_inter(z),Delta_C_Bat_inter(z),C_Rest_V_inter(z)] = Batterie(PWM_inter(z),eta_PWM,I_mot_inter(z),n_Prop,C_Bat,P_Bat_Peukert,Delta_C_Bat_inter(z),t_Flug);
            
            P_Untersuchung(z) = I_Bat_inter(z) * U_Bat_nom;
            
            
            
            
            
            %% Gesamtwirkungsgrad
            
            % Berechnung der induzierten Geschwindigkeiten nach van der Wall
            % (Grundlagen der Hubschrauber-Aerodynamik) (2015) S. 153
            
            vi0 = sqrt(m*g / ( 2*rho(q)*F*n_Prop ) );                          % induzierte Geschwindigkeit im Schwebeflug v_i0
            v = vi0;
            mu_z = -V_A*sin(alpha_inter(z));                                          % Geschwindigkeit durch die Rotorebene
            mu = V_A*cos(alpha_inter(z));                                             % Geschwindigkeit entlang Rotorebene
            krit = 1;
            while krit > 0.0005
                f = v - mu_z - vi0^2 / sqrt(mu^2 + v^2);
                fs = 1 + v * vi0^2 / (mu^2 + v^2)^(3/2);
                v_i_neu = v - f/fs;
                krit = abs(v_i_neu - v) / v_i_neu;
                v = v_i_neu;
            end
            vi_vi0 = (v - mu_z) / vi0;
            vi = vi0 * vi_vi0;                                                  % induzierte Geschwindigkeit im stationaeren Steigflug
            
            % Figure of Merit des Rotors, Bezug auf van der Wall (Grundlagen der Hubschrauber-Aerodynamik) (2015) (S.122)
            eta_ges_inter(z) = (n_Prop * Thrust_inter(z) * (mu_z + vi))/(I_Bat_inter(z) * U_Bat_nom);         % Leistung, die in Schub umgesetzt wird im Verh�ltnis zur aufgebrachten Leistung
            
            
        end
        
        
        
        
        % Werden Grenzen ueberschritten?
        
        % Flugbereichsgrenzen f�r das Fl�chenflugzeug innerhalb der
        % Flugenvellope
        
        if Abfrage_Flugsystem == 0
            
            % aerodynamische Grenze
            % herausnahmen, Annahme Flug mit V* und nicht V_min, i.e. epsilon_opt
            n_z = cos(atan(epsilon));
            V_min = sqrt(2*n_z*m*g/(c_A_plane_max*S*rho(q)));
            
            % Leistungsgrenze
            
            % c_A = C_A0 + alpha * C_Aalpha;
            % c_W = c_W0 + k * C_A^2
            W = V_A^2 * rho(q)/2 * S * c_W_plane;
            T_erf = W;
            
            % Temperaturgrenze
            
            T_zul = 273.15 + t_zul;
            % T_max = 2 * T_zul / ((kappa - 1) * (V_A/a)^2 +2);
            V_max_T = sqrt((T_zul-T)*2/(T*(kappa-1))) * a;
            
            % Begrenzung durch Festigkeit
            
            q_max = rho(q)/2 * V_A^2;
            V_max_q = sqrt(q_zul * 2 /rho(q));
            
            if  V_A > V_max_q || V_A >= V_max_T || V_min > V_A || T_erf > max(max(T_map)) % || T_max > T_zul
                C_Rest_V(q) = NaN;
                Omega(q) = NaN;
                U_mot(q) = NaN;
                I_mot(q) = NaN;
                I_Bat(q) = NaN;
                PWM(q) = NaN;
                eta_ges(q) = NaN;
            end
        end
        
        
        % Wenn Grenzen ueberschritten werden, Resultate entfernen
        
        if C_Rest_V_inter(z) < 0.0 || U_mot_inter(z) > U_Bat_nom || U_mot_inter(z) <= 0 || C_Rate_inter(z) > C_Rate_max || I_mot_inter(z) > I_max || alpha_inter(z) > alpha_stall || M_tip_inter(z) >= 1
            C_Rest_V_inter(z) = NaN;
            Omega_inter(z) = NaN;
            U_mot_inter(z) = NaN;
            I_mot_inter(z) = NaN;
            I_Bat_inter(z) = NaN;
            PWM_inter(z) = NaN;
            eta_ges_inter(z) = NaN;
        end
        
        z = z+1;
        
        
    end
    
    % Kriterium zur Auswahl der optimalen Steiggeschwindigkeit

    if mean(isnan(Delta_C_Bat_inter)) ~= 1                                 % <-- hier noch anderes Krit zur Auswahl aussuchen
        
        % Kriterium f�r optimale Steiggeschwindigkeit
        ind_opt = find(Delta_C_Bat_inter == min(Delta_C_Bat_inter));
%         ind_opt = find(P_Untersuchung == min(P_Untersuchung));			% Index der opt. Steiggeschw. finden
        
        if length(ind_opt) > 1
            ind_opt = ind_opt(1);
        end
        

        % �bergabe der Leistungswerte f�r die optimale Geschwindigkeit und
        % Festlegen f�r entsprechenden H�henschritt
        Thrust(q) = Thrust_inter(ind_opt);
        Theta(q) = Theta_inter(ind_opt);
        alpha(q) = alpha_inter(ind_opt);
        Omega(q) = Omega_inter(ind_opt);
        I_mot(q) = I_mot_inter(ind_opt);
        U_mot(q) = U_mot_inter(ind_opt);
        C_Rate(q) = C_Rate_inter(ind_opt);
        C_Rest_V(q) = C_Rest_V_inter(ind_opt);
        I_Bat(q) = I_Bat_inter(ind_opt);
        tau(q) = tau_inter(ind_opt);
        M_tip(q) = M_tip_inter(ind_opt);
        PWM(q) = PWM_inter(ind_opt);
        eta_ges(q) = eta_ges_inter(ind_opt);
        Delta_C_Bat(q) = Delta_C_Bat_inter(ind_opt);
        
        
        V_Kg_opt(q) = V_Kg(ind_opt);		% Bestimmung opt. Stgeschw. und speichern in Vektor f�r jeden H�henabschnitt
        
    else
        
        V_Kg_opt(q) = NaN;
        Thrust(q) = NaN;
        Theta(q) = NaN;
        alpha(q) = NaN;
        Omega(q) = NaN;
        I_mot(q) = NaN;
        U_mot(q) = NaN;
        C_Rate(q) = NaN;
        C_Rest_V(q) = NaN;
        I_Bat(q) = NaN;
        tau(q) = NaN;
        M_tip(q) = NaN;
        PWM(q) = NaN;
        eta_ges(q) = NaN;
        
        
        
    end 
    
    
    
    H(q) = H_oben;			% Speichern der H�he im Vektor
    q = q+1;                % Erh�hen des Z�hlers um 1
    
    
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot

% Darstellung der Ergenisse in Diagrammen

% Restladung �ber der H�he
figure(figure_C_Rest_V)
plot(H,C_Rest_V*100,'LineWidth',2);
grid on
hold on
xlabel('H�he [m]')
ylabel('Restladung der Batterie [%]')


% Drehzahl �ber der H�he
figure(figure_omega)
plot(H,Omega/(2*pi)*60,'LineWidth',2)
grid on
hold on
xlabel('H�he [m]')
ylabel('RPM')

% Motorstrom �ber der H�he
figure(figure_I_mot)
plot(H,I_mot,'LineWidth',2)
grid on
hold on
xlabel('H�he [m]')
ylabel('I_{mot} [A]')

% Motorspannung �ber der H�he
figure(figure_U_mot)
plot(H,U_mot,'LineWidth',2)
grid on
hold on
xlabel('H�he [m]')
ylabel('U_{mot} [V]')

% Batteriestrom �ber der H�he
figure(figure_I_Bat)
plot(H,I_Bat,'LineWidth',2)
grid on
hold on
xlabel('H�he [m]')
ylabel('I_{Bat} [A]')

% PWM �ber der H�he
figure(figure_PWM)
plot(H,PWM*100,'LineWidth',2)
grid on
hold on
xlabel('H�he [m]')
ylabel('PWM [%]')

% Wirkungsgrad
figure(figure_eta)
plot(H,eta_ges*100,'LineWidth',2)
grid on
hold on
xlabel('H�he [m]')
ylabel('eta_{ges} [%]')

% Steiggeschwindigkeit
figure(figure_vkg)
plot(H,V_Kg_opt,'LineWidth',2)
grid on
hold on
xlabel('H�he [m]')
ylabel('Bahngeschwindigkeit [m/s]')

%% Datei abspeichern
% ImageSizeX = 40;
% ImageSizeY = 30;
% figure(figure_C_Rest_V)
% set(gcf,'PaperUnits','centimeters', 'PaperPosition', [0 0 ImageSizeX ImageSizeY]); 
% set(gcf,'Units','centimeters', 'PaperSize', [ImageSizeX ImageSizeY]); 
% saveas(gcf,'C_Rest_V', 'jpg'); 
% figure(figure_omega)
% saveas(gcf,'omega', 'jpg');
% figure(figure_I_mot)
% saveas(gcf,'I_mot', 'jpg');
% figure(figure_U_mot)
% saveas(gcf,'U_mot', 'jpg');
% figure(figure_I_Bat)
% saveas(gcf,'I_Bat', 'jpg');
% figure(figure_PWM)
% saveas(gcf,'PWM', 'jpg');