%% Leistungsberechnung f�r Flugsysteme %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% *************************************************************************


%% Initialisierungen %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Batterie

U_Bat_nom = N_Bat_cell * U_Bat_cell;        % nominale Batteriespannung
U_Bat_min = N_Bat_cell * U_Bat_cell_min;    % minimale Batteriespannung
% C_Bat = E_Dichte * m_Bat / U_Bat_nom;       % Kapazitaet der Batterie in As
C_Bat = N_Bat_cell_p*C_Bat_cell*3600;
% Delta_C_Bat = 0;                            % Initialisierung Batteriekapazit�t, die nach jedem delta_h gebraucht wird

% Normzelle erzeugen
DATA = Elektromodellflug;

% L�schen der Ausrei�er
DATA(63,:) = [];       % id_bat = 63
DATA(40,:) = [];       % id_bat = 40
DATA(30,:) = [];       % id_bat = 30
DATA(14,:) = [];       % id_bat = 14
DATA(38,:) = [];       % id_bat = 38

[Q,Qnom,Qexp,Vfull,Vexp,Vnom,i] = Normcell(DATA);             % Normzelle generieren
Batterie_data = [Q Qnom Qexp Vfull Vexp Vnom i 0];            % Zwischenbelegung
Cnom = C_Bat/1000;                                            % Nominelle Kapazit�t



%% Propeller

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





%% Umgebungsparameter

T_11 = T_0 - 0.0065 *(11000-H_0);                   % T in 11000m H�he
rho_11 = rho_0 * (1 - 0.0065*(11000/T_0))^4.256;    % Dichte in 11000m H�he
p_11 = p_0 * (1 - 0.0065*(11000/T_0))^5.256;        % Druck in 11000m H�he



%% Intialisierung der Matrizen f�r jeden H�henabschnitt

lengthi = floor(abs(H_max - H_0) / Delta_H + 1);
lengthgamma = floor(abs(gamma_max - gamma_min) / gamma_Delta + 1);

H = zeros(lengthi,1);
U_mot = zeros(lengthi,1);
I_mot = zeros(lengthi,1);
Thrust = zeros(lengthi,1);
Omega = zeros(lengthi,1);
tau = zeros(lengthi,1);
C_Rate = zeros(lengthi,1);
alpha = zeros(lengthi,1);
C_Rest_V = zeros(lengthi,1);
Theta = zeros(lengthi,1);
rho = zeros(lengthi,1);
I_Bat = zeros(lengthi,1);
% U_Bat = zeros(lengthi,1);
P_Bat = zeros(lengthi,1);
PWM = zeros(lengthi,1);
M_tip = zeros(lengthi,1);
eta_prop = zeros(lengthi,1);
eta_ges = zeros(lengthi,1);
gamma_Flaechenflzg = zeros(lengthi,1);
Flugzustand_Flaechenflzg = zeros(lengthi,1);
Delta_C_Bat = zeros(lengthi,1);


% U_Bat(1) = U_Bat_nom;
i_int = zeros(lengthi+1,1);
i_int_inter = zeros(lengthgamma,1);





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Programmanfang

x = 1;                                      % Z�hler intialisiern
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
    
    rho(x) = rho_unten + (rho_oben - rho_unten)/2;                  % Berechnung der mittleren Dichte im Intervall
    p = (p_oben + p_unten)/2;                                       % mittlerer Druck im Intervall
    a = sqrt(kappa*p/rho(x));                                       % Schallgeschwindigkeit in m/s
    
    
    % Dichte an der oberen (_2) und unteren (_1) Intervallgrenze
    if x == 1
        rho_1 = rho_0;
        rho_2 = rho(x);
    else
        rho_1 = rho(x-1);
        rho_2 = rho(x);
    end
    
    
    T_map = T_map * rho(x)/rho_1;                                   % Anpassung des Schubkennfeldes an die sich �ndernde Dichte
    P_map = P_map * rho(x)/rho_1;                                   % Anpassung des Leistungskennfeldes an die sich �ndernde Dichte
    TAU_map = TAU_map * rho(x)/rho_1;                               % Anpassung des Drehmomentkennfeldes an die sich �ndernde Dichte
    
    %% Beginn der Leistungsberechnung f�r Flugsysteme
    
    % Initialisierung des Auswahlvektors f�r den Steigflug
    
    
    Thrust_inter = zeros(lengthgamma, 1);
    Theta_inter = zeros(lengthgamma, 1);
    alpha_inter = zeros(lengthgamma, 1);
    Omega_inter = zeros(lengthgamma, 1);
    I_mot_inter = zeros(lengthgamma, 1);
    U_mot_inter = zeros(lengthgamma, 1);
    C_Rate_inter = zeros(lengthgamma, 1);
    C_Rest_V_inter = zeros(lengthgamma, 1);
    I_Bat_inter = zeros(lengthgamma, 1);
    U_bat_inter = zeros(lengthgamma, 1);
    tau_inter = zeros(lengthgamma, 1);
    M_tip_inter = zeros(lengthgamma, 1);
    PWM_inter = zeros(lengthgamma, 1);
    eta_ges_inter = zeros(lengthgamma, 1);
    Flugzustand_Flaechenflzg_inter = zeros(lengthgamma,1);
    gamma_Flaechenflzg_inter = zeros(lengthgamma,1);
    Delta_C_Bat_inter = zeros(lengthgamma,1);
    
    Bestimmung_gamma = zeros(lengthgamma, 1);
    
    
    z = 1;
    
    for gamma_variabel = gamma_min:gamma_Delta:gamma_max	    % Variation des Bahnneigungswinkel f�r das Fl�chenflugzeug
        
        gamma_Flaechenflzg_inter(z) = gamma_variabel;
        
        if Abfrage_Flugsystem == 1
            
            % MULTICOPTER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            m = m_copter + m_Bat + m_Mot * n_Prop + m_nutz;                     % Gesamtmasse des Quadrocopters
            
            t_Flug = Delta_H / V_Kg;                                            % Flugzeit
            
            
            % Aerodynamik berechnen
            
            [Thrust_inter(z),Theta_inter(z),V_A,alpha_inter(z)] = MulticopterAerodynamik(u_Wg,V_Kg,gamma_copter,c_W_copter_seitlich,c_W_copter_oben,c_A_copter_max,rho(x),A_copter,m,g);
            
            
        else
            
            
            % FL�CHENFLUGZEUG %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            m = m_flugzeug + m_Bat + m_Mot * n_Prop + m_nutz;                   % Gesamtmasse des Fl�chenflugzeugs
            
            % Aerodynamik
            
            [Thrust_inter(z),V_A,Flugzustand_Flaechenflzg_inter(z)] = FlaechenflugzeugAerodynamik(m,g,E_stern,V_stern,rho_stern,E,gamma_variabel,rho(x));
            
            V_H = V_A * sind(gamma_variabel);
            
            t_Flug = Delta_H / V_H;
            
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
            
            [PWM_inter(z),eta_PWM] = ESC(U_mot_inter(z),U_Bat_nom);			% <-- hier U_bat_inter
            
            
            % Batteriezustand berechnen
            
            Delta_C_Bat_inter(z) = Delta_C_Bat(x);                           % Anpassung der Batteriekapazit�t
            [I_Bat_inter(z),C_Rate_inter(z),Delta_C_Bat_inter(z),C_Rest_V_inter(z)] = Batterie(PWM_inter(z),eta_PWM,I_mot_inter(z),n_Prop,C_Bat,P_Bat_Peukert,Delta_C_Bat_inter(z),t_Flug);
            
            %		i_int_inter(z) = i_int(x);
            
            % 		[I_Bat(x),U_bat(x),C_rate(x),Delta_C_bat,C_Rest_V(x),i_int_inter(z)] = Batterie(Batterie_data,Cnom,PWM(x),...
            %           		eta_PWM,n_Prop,i_int_inter(z),U_bat(x),C_bat,Delta_C_bat,I_mot(x),N_bat_cell,P_Bat_Peukert,t_Flug)
            
            
            % Gesamtwirkungsgrad
            
            % Berechnung der induzierten Geschwindigkeiten nach van der Wall
            % (Grundlagen der Hubschrauber-Aerodynamik) (2015) S. 153
            
            vi0 = sqrt(m*g / ( 2*rho(x)*F*n_Prop ) );                          % induzierte Geschwindigkeit im Schwebeflug v_i0
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
            %        	eta_prop(x) = (Thrust(x) * (V_A + vi))/(tau(x) .* Omega(x));
            
            eta_ges_inter(z) = (n_Prop * Thrust_inter(z) * (mu_z + vi))/(I_Bat_inter(z) * U_Bat_nom);         % Leistung, die in Schub umgesetzt wird im Verh�ltnis zur aufgebrachten Leistung
            
            
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
        
        
        
        if Abfrage_Flugsystem == 1
            
            break;								% Abbruch f�r den Fall eines Multicopters um for-Schleife zu verlassen
            
        else
            
            Bestimmung_gamma(z) = Delta_C_Bat_inter(z) * U_Bat_nom;		% Berechnung der aufgebrachten Energiemenge
            
        end
        
        
        z = z+1;			% Erh�hung der Z�hlervariablen f�r die gamma-Schleife
        
    end
    
    
    % Hier Auswahl des Steigwinkels f�r Fl�chenflugzeug
    
    
    if Abfrage_Flugsystem == 1
        
        Thrust(x) = Thrust_inter(1);
        Theta(x) = Theta_inter(1);
        alpha(x) = alpha_inter(1);
        Omega(x) = Omega_inter(1);
        I_mot(x) = I_mot_inter(1);
        U_mot(x) = U_mot_inter(1);
        C_Rate(x) = C_Rate_inter(1);
        C_Rest_V(x) = C_Rest_V_inter(1);
        I_Bat(x) = I_Bat_inter(1);
        tau(x) = tau_inter(1);
        M_tip(x) = M_tip_inter(1);
        PWM(x) = PWM_inter(1);
        eta_ges(x) = eta_ges_inter(1);
        Delta_C_Bat(x) = Delta_C_Bat(x) + Delta_C_Bat_inter(1);
        % 	i_int(x) = i_int_inter(1);
        
    else
        
        if mean(isnan(Delta_C_Bat_inter)) ~= 1                                 % <-- hier noch anderes Krit zur Auswahl aussuchen
            
            % Kriterium f�r optimalen Steigwinkel
            ind_opt = find(Bestimmung_gamma == min(Bestimmung_gamma));
            
            
            if length(ind_opt) > 1
                ind_opt = ind_opt(1);
            end
            
            
            % �bergabe der Leistungswerte f�r die optimale Geschwindigkeit und
            % Festlegen f�r entsprechenden H�henschritt
            Thrust(x) = Thrust_inter(ind_opt);
            Theta(x) = Theta_inter(ind_opt);
            alpha(x) = alpha_inter(ind_opt);
            Omega(x) = Omega_inter(ind_opt);
            I_mot(x) = I_mot_inter(ind_opt);
            U_mot(x) = U_mot_inter(ind_opt);
            C_Rate(x) = C_Rate_inter(ind_opt);
            C_Rest_V(x) = C_Rest_V_inter(ind_opt);
            I_Bat(x) = I_Bat_inter(ind_opt);
            tau(x) = tau_inter(ind_opt);
            M_tip(x) = M_tip_inter(ind_opt);
            PWM(x) = PWM_inter(ind_opt);
            eta_ges(x) = eta_ges_inter(ind_opt);
            Delta_C_Bat(x) = Delta_C_Bat(x) + Delta_C_Bat_inter(ind_opt);   % �nderung
            i_int(x) = i_int_inter(ind_opt);
            Flugzustand_Flaechenflzg(x) = Flugzustand_Flaechenflzg_inter(ind_opt);
            
            
            gamma_Flaechenflzg(x) = gamma_Flaechenflzg_inter(ind_opt);		% Bestimmung opt. Stgeschw. und speichern in Vektor f�r jeden H�henabschnitt
            
        else
            
            Thrust(x) = NaN;
            Theta(x) = NaN;
            alpha(x) = NaN;
            Omega(x) = NaN;
            I_mot(x) = NaN;
            U_mot(x) = NaN;
            C_Rate(x) = NaN;
            C_Rest_V(x) = NaN;
            I_Bat(x) = NaN;
            tau(x) = NaN;
            M_tip(x) = NaN;
            PWM(x) = NaN;
            eta_ges(x) = NaN;
            
            
            
        end
        
    end
    
    
    H(x) = H_oben;			% Speichern der H�he im Vektor
    x = x+1;				% Erh�hung der Z�hlervariablen f�r die H�hen-Schleife
    
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

% Steigwinkel Flaechenflugzeug
figure(figure_gamma)
for i = 1:length(gamma_Flaechenflzg)
    
    if Flugzustand_Flaechenflzg == 0
        
        plot(gamma_Flaechenflzg(i),H(i),'b.');
        grid on
        hold on
        
    elseif Flugzustand_Flaechenflzg == 1
        
        plot(gamma_Flaechenflzg(i),H(i),'y.');
        grid on
        hold on
        
    else
        
        plot(gamma_Flaechenflzg(i),H(i),'r.');
        grid on
        hold on
        
    end
    
end
xlabel('Bahnneigungswinkel [�]')
ylabel('H�he [m]')



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