%% Leistungsberechnung für Flugsysteme %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% *************************************************************************


%% Initialisierungen %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Batterie

U_Bat_nom = N_Bat_cell * U_Bat_cell;        % nominale Batteriespannung
U_Bat_min = N_Bat_cell * U_Bat_cell_min;    % minimale Batteriespannung
C_Bat = E_Dichte * m_Bat / U_Bat_nom;       % Kapazitaet der Batterie in As
% C_Bat = N_Bat_cell_p*C_Bat_cell*3600;       % Kapazität in As
% Delta_C_Bat = 0;                            % Initialisierung Batteriekapazität, die nach jedem delta_h gebraucht wird

% Normzelle erzeugen
DATA = Elektromodellflug;

% Löschen der Ausreißer
DATA(63,:) = [];       % id_bat = 63
DATA(40,:) = [];       % id_bat = 40
DATA(30,:) = [];       % id_bat = 30
DATA(14,:) = [];       % id_bat = 14
DATA(38,:) = [];       % id_bat = 38

[Q,Qnom,Qexp,Vfull,Vexp,Vnom,i] = Normcell(DATA);             % Normzelle generieren
Batterie_data = [Q Qnom Qexp Vfull Vexp Vnom i 0];            % Zwischenbelegung
Cnom = C_Bat/1000;                                            % Nominelle Kapazität



%% Propeller

%Entnahme des Durchmessers und des Pitches aus dem Propellernahmen
D = str2double(prop_name(1:strfind(prop_name,'x')-1));      % Durchmesser extrahieren
P_75 = prop_name(strfind(prop_name,'x')+1:end);             % Pitch extrahieren
while isnan(str2double(P_75)) == 1
    P_75(end) = [];
end
P_75 = str2double(P_75);                    % Pitch festlegen
R = D * 0.0254 / 2;                         % Propellerradius in Meter
F = pi * R^2;                               % Fläche eines Propellers in Quadratmeter
Theta_75 = atan( 4*P_75 / (3*pi * D) );     % geometrischer Anstellwinkel des Propellers bei 75% des Radius
[RPM_map, V_map, T_map, P_map, TAU_map] = Propeller_map(DATA_APC,prop_name);    % Aufbau des Kennfeldes





%% Umgebungsparameter

T_11 = T_0 - 0.0065 *(11000-H_0);                   % T in 11000m Höhe
rho_11 = rho_0 * (1 - 0.0065*(11000/T_0))^4.256;    % Dichte in 11000m Höhe
p_11 = p_0 * (1 - 0.0065*(11000/T_0))^5.256;        % Druck in 11000m Höhe



%% Intialisierung der Matrizen für jeden Höhenabschnitt

lengthi = floor(abs(H_max - H_0) / Delta_H + 1);
lengthgamma = floor(abs(gamma_max - gamma_min) / gamma_Delta + 1);
% v_ver_max = ceil(V_stern);            
% lengthb_vert = floor(abs(v_ver_max - v_vert_min) / v_vert_Delta + 1);

% Umgebung
H = zeros(lengthi,1);
rho = zeros(lengthi,1);
% Multicopter
alpha = zeros(lengthi,1);
Theta = zeros(lengthi,1);
% Flächenflugzeug
gamma_Flaechenflzg = zeros(lengthi,1);
Flugzustand_Flaechenflzg = zeros(lengthi,1);
V_Flaechenflugzeug = zeros(lengthi,1);
% Propeller
Thrust = zeros(lengthi,1);
Omega = zeros(lengthi,1);
tau = zeros(lengthi,1);
M_tip = zeros(lengthi,1);
% Motor
U_mot = zeros(lengthi,1);
I_mot = zeros(lengthi,1);
% ESC
PWM = zeros(lengthi,1);
% Batterie
C_Rate = zeros(lengthi,1);
C_Rest_V = zeros(lengthi,1);
I_Bat = zeros(lengthi,1);
U_Bat = zeros(lengthi+1,1);
U_Bat(1) = U_Bat_nom;                               % Startspannung ist die Nominalspannung
P_Bat = zeros(lengthi,1);
Delta_C_Bat = zeros(lengthi+1,1);
i_int = zeros(lengthi+1,1);
% Gesamtsystem
eta_prop = zeros(lengthi,1);
eta_ges = zeros(lengthi,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Programmanfang

x = 1;                                      % Zähler intialisiern
for h_variabel = H_0:Delta_H:H_max
    
    %% Umgebungsparameter als Funktion der Höhe
    % Berechnung der Flughöhe für iterative Schritte und der mittleren
    % Dichte zwischen den Diskretisierungspunkten
    
    H_unten = h_variabel;
    H_oben = h_variabel + Delta_H;
    H_mitte = (H_oben + H_unten)/2;
    
    % Berechnung der Dichte an den Intevallgrenzen nach Normatmosphärenbedingungen
    
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
    
    
    T_map = T_map * rho(x)/rho_1;                                   % Anpassung des Schubkennfeldes an die sich ändernde Dichte
    P_map = P_map * rho(x)/rho_1;                                   % Anpassung des Leistungskennfeldes an die sich ändernde Dichte
    TAU_map = TAU_map * rho(x)/rho_1;                               % Anpassung des Drehmomentkennfeldes an die sich ändernde Dichte
    
    
    %% Beginn der Leistungsberechnung für Flugsysteme
    
    % Initialisierung des Auswahlvektors für den Steigflug
 
    % Gesamtlänge der Vektoren: lengthall = lengthgamma + lengthb_vert;
    
    
    % Multicopter
    Theta_inter = zeros(lengthgamma, 1);
    alpha_inter = zeros(lengthgamma, 1);
    % Flaechenflugzeug
    Flugzustand_Flaechenflzg_inter = zeros(lengthgamma,1);
    gamma_Flaechenflzg_inter = zeros(lengthgamma,1);
    Bestimmung_gamma = zeros(lengthgamma, 1);
    V_Flaechenflugzeug_inter = zeros(lengthgamma, 1);
    % Propeller
    Omega_inter = zeros(lengthgamma, 1);
    tau_inter = zeros(lengthgamma, 1);
    M_tip_inter = zeros(lengthgamma, 1);
    % Motor
    I_mot_inter = zeros(lengthgamma, 1);
    U_mot_inter = zeros(lengthgamma, 1);
    % ESC
    PWM_inter = zeros(lengthgamma, 1);
    % Batterie
    I_Bat_inter = zeros(lengthgamma, 1);
    U_Bat_inter = zeros(lengthgamma, 1);
    C_Rate_inter = zeros(lengthgamma, 1);
    C_Rest_V_inter = zeros(lengthgamma, 1);
    Delta_C_Bat_inter = zeros(lengthgamma,1);
    i_int_inter = zeros(lengthgamma,1);
    % Gesamtsystem
    Thrust_inter = zeros(lengthgamma, 1);
    eta_ges_inter = zeros(lengthgamma, 1);
    
%   b = 0;                          % Kontrollfaktor
%   V_vert = 0;                     % Vertikalgeschwindigkeit
    
    z = 1;
    
    for gamma_variabel = gamma_min:gamma_Delta:gamma_max	    % Variation des Bahnneigungswinkel für das Flächenflugzeug
        
        gamma_Flaechenflzg_inter(z) = gamma_variabel;
        
        if Abfrage_Flugsystem == 1
            
            % MULTICOPTER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            m = m_copter + m_Bat + m_Mot * n_Prop + m_nutz;                     % Gesamtmasse des Quadrocopters
            t_Flug = Delta_H / V_Kg;                                            % Flugzeit
            
            % Aerodynamik
            [Thrust_inter(z),Theta_inter(z),V_A,alpha_inter(z)] = MulticopterAerodynamik(u_Wg,V_Kg,gamma_copter,c_W_copter_seitlich,c_W_copter_oben,c_A_copter_max,rho(x),A_copter,m,g);
            
            
        else
            
            
            % FLÄCHENFLUGZEUG %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            m = m_flugzeug + m_Bat + m_Mot * n_Prop + m_nutz;                   % Gesamtmasse des Flächenflugzeugs
            
            % Hier fehlt noch Behandlung des Zustandes, wenn kein
            % Flugzustand == 2 erreicht wird
            
            
%             % Aerodynamik
%              if b == 0             
%                   [Thrust_inter(z),V_A,Flugzustand_Flaechenflzg_inter(z)] = FlaechenflugzeugAerodynamik(m,g,E_stern,V_stern,rho_stern,E,gamma_variabel,rho(x));  
%                   V_H = V_A * sind(gamma_variabel);
%                   t_Flug = Delta_H / V_H;
%                   V_Flaechenflugzeug_inter(z) = V_A;
%              end


%                 if Flugzustand_Flaechenflzg == 2 || b == 1 || gamma_variabel >= 90
%                     b = 1;
%                     V_vert = V_vert + v_Delta;
%                     gamma_Flaechenflzg_inter(z) = 90;
%                     [Thrust_inter(z)] = Vertikalflug(m,g,E_stern,V_stern,rho_stern,V,rho(x));
%                     t_Flug = Delta_H / V_vert;
%                     V_Flaechenflzg_inter(z) = V_vert;
%                 elseif b == 0 && gamma_variabel >= 90
%                     Thrust = 0;
                      
%                 end     
%

%             alpha_inter(z) = - 90 * pi/180;
        end
        
        
        Thrust_inter(z) = Thrust_inter(z) / n_Prop;                        % Schub auf n Propeller verteilen
        
        

        if Thrust_inter(z) > max(max(T_map))                               % wenn Schub zu gross (Ergebnis verwerfen)
            
            Omega_inter(z) = NaN;
            I_mot_inter(z) = NaN;
            C_Rate_inter(z) = NaN;
            C_Rest_V_inter(z) = NaN;
            
            
        else

            % Drehzahl und Drehmoment bestimmen
                    
            [Omega_inter(z),tau_inter(z)] = Propeller(V_A, alpha_inter(z), Thrust_inter(z), RPM_map, V_map, T_map, TAU_map);
            
            
            % Wie groß ist die Blattspitzengeschwindigkeit?
            M_tip_inter(z) = (Omega_inter(z) * R)/a;                       % Blattspitzengeschwindigkeit in Ma
            
            
            % Motorzustand berechnen
            [U_mot_inter(z),I_mot_inter(z)] = Motor(tau_inter(z),K_V,I_0,R_i,Omega_inter(z));
            
            
            % Zustand der Motorregler berechnen
            [PWM_inter(z),eta_PWM] = ESC(U_mot_inter(z),U_Bat(x));         % <-- hier U_bat_inter
            
            
            % Batteriezustand berechnen
            Delta_C_Bat_inter(z) = Delta_C_Bat(x);                         % Anpassung der Batteriekapazität
            i_int_inter(z) = i_int(x);                                     % Übergabe des Integrals der Spannung vom letzten Schritt
            %            [I_Bat_inter(z),C_Rate_inter(z),Delta_C_Bat_inter(z),C_Rest_V_inter(z)] = Batterie(PWM_inter(z),eta_PWM,I_mot_inter(z),n_Prop,C_Bat,P_Bat_Peukert,Delta_C_Bat_inter(z),t_Flug);
            
            [I_Bat_inter(z),U_Bat_inter(z),C_Rate_inter(z),Delta_C_Bat_inter(z),C_Rest_V_inter(z),i_int_inter(z)] = Batterie_neu(Batterie_data,...
                Cnom,PWM_inter(z),eta_PWM,n_Prop,i_int_inter(z),U_Bat_inter(z),C_Bat,Delta_C_Bat_inter(z),I_mot_inter(z),N_Bat_cell,P_Bat_Peukert,t_Flug);
            
            %% Gesamtwirkungsgrad
            
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
            
            eta_ges_inter(z) = (n_Prop * Thrust_inter(z) * (mu_z + vi))/(I_Bat_inter(z) * U_Bat_nom);         % Leistung, die in Schub umgesetzt wird im Verhältnis zur aufgebrachten Leistung
            
            
        end
        
        % Wenn Grenzen ueberschritten werden, Resultate entfernen
        alpha_inter(z) = alpha_inter(z)*180/pi;
        
        if C_Rest_V_inter(z) < 0.0 || U_mot_inter(z) > U_Bat(x) || U_mot_inter(z) <= 0 || C_Rate_inter(z) > C_Rate_max || I_mot_inter(z) > I_max || ...
                alpha_inter(z) > alpha_stall || M_tip_inter(z) >= 1 || I_Bat_inter(z) <= 0     % ||  PWM_inter(z) > 1.0
            C_Rest_V_inter(z) = NaN;
            Omega_inter(z) = NaN;
            U_mot_inter(z) = NaN;
            I_mot_inter(z) = NaN;
            I_Bat_inter(z) = NaN;
            PWM_inter(z) = NaN;
            eta_ges_inter(z) = NaN;
            gamma_Flaechenflzg_inter(z) = NaN;
            V_Flaechenflugzeug_inter(z) = NaN;
            
        end
        
        % Muss der Steigwinkel variiert werden?
        
        if Abfrage_Flugsystem == 1
            break;								% Abbruch für den Fall eines Multicopters um for-Schleife zu verlassen, da gamma vorgegeben
        else
            Bestimmung_gamma(z) = Delta_C_Bat_inter(z) * U_Bat_inter(z);        % Berechnung der aufgebrachten Energiemenge
        end
        
        
        z = z+1;			% Erhöhung der Zählervariablen für die gamma-Schleife
        
    end

    if Abfrage_Flugsystem == 1
        
        % Multicopter
        Theta(x) = Theta_inter(1);
        alpha(x) = alpha_inter(1);
        % Propeller
        Thrust(x) = Thrust_inter(1);
        Omega(x) = Omega_inter(1);
        tau(x) = tau_inter(1);
        M_tip(x) = M_tip_inter(1);
        % Motor
        I_mot(x) = I_mot_inter(1);
        U_mot(x) = U_mot_inter(1);
        % ESC
        PWM(x) = PWM_inter(1);
        % Batterie
        I_Bat(x) = I_Bat_inter(1);
        U_Bat(x+1) = U_Bat_inter(1);
        C_Rate(x) = C_Rate_inter(1);
        C_Rest_V(x) = C_Rest_V_inter(1);
        Delta_C_Bat(x+1) = Delta_C_Bat_inter(1);
        i_int(x+1) = i_int_inter(1);
        % Gesamtsystem
        eta_ges(x) = eta_ges_inter(1);
        
    else
        if isnan(nanmean(Omega_inter)) ~= 1 && isnan(nanmean(I_mot_inter)) ~= 1 &&  isnan(nanmean(U_mot_inter)) ~= 1 && ...     % Wenn alle der Vektoren nicht nur NaN
                isnan(nanmean(PWM_inter)) ~= 1 && isnan(nanmean(C_Rest_V_inter)) ~= 1 && isnan(nanmean(I_Bat_inter)) ~= 1       % enthalten (unfliegbarer Zustand) 
            
            % Kriterium für optimalen Steigwinkel
            y = 1;
            while y < length(Bestimmung_gamma)                             % Suche und finde optimalen Steigwinkel
                A = Bestimmung_gamma;
                A = sort(A);
                ind_opt = find(Bestimmung_gamma == A(y));
                
                if length(ind_opt) > 1
                    ind_opt = ind_opt(1);
                end
                
                if ~isnan(Thrust_inter(ind_opt)) && ~isnan(Omega_inter(ind_opt)) && ~isnan(I_mot_inter(ind_opt)) && ~isnan(U_mot_inter(ind_opt)) && ~isnan(PWM_inter(ind_opt)) ...
                        && ~isnan(I_Bat_inter(ind_opt)) && ~isnan(U_Bat_inter(ind_opt)) && ~isnan(C_Rate_inter(ind_opt)) && ~isnan(C_Rest_V_inter(ind_opt))
                    
                    break;                                                 % Unterbreche Schleife, falls alle Leistungsparameter physikalisch realistische Werte besitzen
                end
                
                y = y+1;
            end
            
            % Übergabe der Leistungswerte für die optimale Geschwindigkeit und
            % Festlegen für entsprechenden Höhenschritt
            
            % Flächenflugzeug
            Theta(x) = Theta_inter(ind_opt);
            alpha(x) = alpha_inter(ind_opt);
            gamma_Flaechenflzg(x) = gamma_Flaechenflzg_inter(ind_opt);	    % Bestimmung opt. Stgeschw. und speichern in Vektor für jeden Höhenabschnitt
            Flugzustand_Flaechenflzg(x) = Flugzustand_Flaechenflzg_inter(ind_opt);
            V_Flaechenflugzeug(x) = V_Flaechenflugzeug_inter(ind_opt);
            % Propeller
            Thrust(x) = Thrust_inter(ind_opt);
            Omega(x) = Omega_inter(ind_opt);
            tau(x) = tau_inter(ind_opt);
            M_tip(x) = M_tip_inter(ind_opt);
            % Motor
            I_mot(x) = I_mot_inter(ind_opt);
            U_mot(x) = U_mot_inter(ind_opt);
            % ESC
            PWM(x) = PWM_inter(ind_opt);
            % Batterie
            I_Bat(x) = I_Bat_inter(ind_opt);
            U_Bat(x+1) = U_Bat_inter(ind_opt);
            C_Rate(x) = C_Rate_inter(ind_opt);
            C_Rest_V(x) = C_Rest_V_inter(ind_opt);
            Delta_C_Bat(x+1) = Delta_C_Bat_inter(ind_opt);
            i_int(x+1) = i_int_inter(ind_opt);
            % Gesamtsystem
            eta_ges(x) = eta_ges_inter(ind_opt);
            
        else
            
            % Ansonsten verwerfe alle Werte
            % Flächenflugzeug
            Theta(x) = NaN;
            alpha(x) = NaN;
            gamma_Flaechenflzg(x) = NaN;
            V_Flaechenflugzeug(x) = NaN;
            % Propeller
            Thrust(x) = NaN;
            Omega(x) = NaN;
            tau(x) = NaN;
            M_tip(x) = NaN;
            % Motor
            I_mot(x) = NaN;
            U_mot(x) = NaN;
            % ESC
            PWM(x) = NaN;
            % Batterie
            I_Bat(x) = NaN;
            U_Bat(x+1) = NaN;
            C_Rate(x) = NaN;
            C_Rest_V(x) = NaN;
            % Gesamtsystem
            eta_ges(x) = NaN;
            
        end
        
    end
    
    
    H(x) = H_oben;			% Speichern der Höhe im Vektor
    x = x+1;				% Erhöhung der Zählervariablen für die Höhen-Schleife
    
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot

% Darstellung der Ergenisse in Diagrammen
figure(figure_ges)
if Abfrage_Flugsystem == 1
    % MULICOPTER
    subplot(421), plot(H,C_Rest_V*100,'LineWidth',2), grid, title('Restladung'), xlabel('Höhe [m]'),ylabel('C_{Bat,Rest} [%]')
    subplot(422), plot(H,Omega/(2*pi)*60,'LineWidth',2), grid, title('Drehzahl'), xlabel('Höhe [m]'),ylabel('Drehzahl [RPM]')
    subplot(423), plot(H,I_mot,'LineWidth',2), grid, title('Motorstrom'), xlabel('Höhe [m]'),ylabel('I_{Mot} [A]')
    subplot(424), plot(H,U_mot,'LineWidth',2), grid,  title('Motorspannung'), xlabel('Höhe [m]'),ylabel('U_{mot} [V]')
    subplot(425), plot(H,I_Bat,'LineWidth',2), grid, title('Batteriestrom'), xlabel('Höhe [m]'),ylabel('I_{Bat} [%]')
    H2 = [0;H];
    subplot(426), plot(H2,U_Bat,'LineWidth',2), grid, title('Batteriespannung'), xlabel('Höhe [m]'),ylabel('U_{Bat} [A]')
    subplot(427), plot(H,PWM*100,'LineWidth',2), grid, title('Pulsweitenmodulation'), xlabel('Höhe [m]'),ylabel('PWM [%]')
    subplot(428), plot(H,eta_ges*100,'LineWidth',2), grid, title('Gesamtwirkungsgrad'), xlabel('Höhe [m]'),ylabel('eta_{ges} [%]')
        
else   
    % FLAECHENFLUGZEUG
    subplot(521), plot(H,C_Rest_V*100,'LineWidth',2), grid, title('Restladung'), xlabel('Höhe [m]'),ylabel('C_{Bat,Rest} [%]')
    subplot(522), plot(H,Omega/(2*pi)*60,'LineWidth',2), grid, title('Drehzahl'), xlabel('Höhe [m]'),ylabel('Drehzahl [RPM]')
    subplot(523), plot(H,I_mot,'LineWidth',2), grid, title('Motorstrom'), xlabel('Höhe [m]'),ylabel('I_{Mot} [A]')
    subplot(524), plot(H,U_mot,'LineWidth',2), grid,  title('Motorspannung'), xlabel('Höhe [m]'),ylabel('U_{mot} [V]')
    subplot(525), plot(H,I_Bat,'LineWidth',2), grid, title('Batteriestrom'), xlabel('Höhe [m]'),ylabel('I_{Bat} [%]')
    H2 = [0;H];
    subplot(526), plot(H2,U_Bat,'LineWidth',2), grid, title('Batteriespannung'), xlabel('Höhe [m]'),ylabel('U_{Bat} [A]')
    subplot(527), plot(H,PWM*100,'LineWidth',2), grid, title('Pulsweitenmodulation'), xlabel('Höhe [m]'),ylabel('PWM [%]')
    subplot(528), plot(H,eta_ges*100,'LineWidth',2), grid, title('Gesamtwirkungsgrad'), xlabel('Höhe [m]'),ylabel('eta_{ges} [%]')
    subplot(529), plot(H,gamma_Flaechenflzg,'LineWidth',2), title('Bahnneigungswinkel'), grid, xlabel('Höhe [m]'),ylabel('gamma [°]')
    subplot(5,2,10), plot(H,V_Flaechenflugzeug,'LineWidth',2), title('absolute Fluggeschwindigkeit'), grid, xlabel('Höhe [m]'),ylabel('V_A [m/s]')
    
end
ImageSizeX = 14;
ImageSizeY = 24;
figure(figure_ges)
set(gcf,'PaperUnits','centimeters', 'PaperPosition', [0 0 ImageSizeX ImageSizeY]); 
set(gcf,'Units','centimeters', 'PaperSize', [ImageSizeX ImageSizeY]); 
saveas(gcf,Dateiname, 'pdf');  

% figure(figure_gamma)
% for i = 1:length(gamma_Flaechenflzg)
% 
%     if Flugzustand_Flaechenflzg == 0
% 
%         plot(H(i),gamma_Flaechenflzg(i),'b.','LineWidth',2);
%         grid on
%         hold on
% 
%     elseif Flugzustand_Flaechenflzg == 1
% 
%         plot(H(i),gamma_Flaechenflzg(i),'y.','LineWidth',2);
%         grid on
%         hold on
% 
%     else
% 
%         plot(H(i),gamma_Flaechenflzg(i),'r.','LineWidth',2);
%         grid on
%         hold on
% 
%     end
% 
% end
% xlabel('Höhe [m]')
% ylabel('Bahnneigungswinkel [°]')


%% Datei abspeichern
% ImageSizeX = 40;
% ImageSizeY = 30;
% set(gcf,'PaperUnits','centimeters', 'PaperPosition', [0 0 ImageSizeX ImageSizeY]);
% set(gcf,'Units','centimeters', 'PaperSize', [ImageSizeX ImageSizeY]);
% figure(figure_C_Rest_V)
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

  %% Datei abspeichern
% ImageSizeX = 14;
% ImageSizeY = 24;
% figure(figure_C_Rest_V)
% figure(figure_omega)
% set(gcf,'PaperUnits','centimeters', 'PaperPosition', [0 0 ImageSizeX ImageSizeY]); 
% set(gcf,'Units','centimeters', 'PaperSize', [ImageSizeX ImageSizeY]); 
% saveas(gcf,Dateiname, 'pdf');  