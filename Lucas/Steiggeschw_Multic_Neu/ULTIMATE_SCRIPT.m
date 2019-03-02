%% Leistungsberechnung für Flugsysteme %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
close all
clear
% *************************************************************************

load('DATA_APC.mat');
load('Elektromodellflug');
load('axi_motor_db.mat');

%% Festlegung des Dateinamen
Dateiname = 'Anz_Prop';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% *************************************************************************
Abfrage_mot = [10 21 28 31];
Abfrage_prop = {'9x4','11x3','11x3','11x3'};
Abfrage_n_prop = [1 2 4 6 8];
Abfrage_c_W = [0.1 0.5 1 1.5 2];
Abfrage_n_bat = [2 4 5 6 8];

% figures definieren
figure_ges = figure;

% Diskretisierung der Steiggeschwindigkeit
V_Kg_min = 1;			% kleinster Bahnneigungswinkel
V_Kg_Delta = 1;		% Schrittweite Batteriemasse
V_Kg_max = 40;			% größter Bahnneigungswinkel

%% Flugparameter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gamma_copter = 90 * pi/180;                 % Bahnanstellwinkel für den Multicopter    

%% Umgebungsparameter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

g = 9.81;                                   % Erdbeschleunigung in m/s^2

H_0 = 0;                                    % Höhe des Abflugplatzes über Normalnull in m
Delta_H = 100;                              % Inkrementweite in m 
H_max = 20000;                              % Maximalhöhe in m

T_0 = 288.15;                               % Temperatur in K am Flugplatz
p_0 = 101325;                               % Druck am Abflugplatz in Pa
rho_0 = 1.225;                              % Dichte am Startort in kg/m^3
kappa = 1.4;                                % Adiabatenexponent

u_Wg = 10;                                  % Seitenwindgeschwindigkeit in m/s

%% Initialisierungen %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Batterie

% U_Bat_nom = N_Bat_cell * U_Bat_cell;        % nominale Batteriespannung
% U_Bat_min = N_Bat_cell * U_Bat_cell_min;    % minimale Batteriespannung
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

%% Umgebungsparameter

T_11 = T_0 - 0.0065 *(11000-H_0);                   % T in 11000m Höhe
rho_11 = rho_0 * (1 - 0.0065*(11000/T_0))^4.256;    % Dichte in 11000m Höhe
p_11 = p_0 * (1 - 0.0065*(11000/T_0))^5.256;        % Druck in 11000m Höhe



%% Intialisierung der Matrizen für jeden Höhenabschnitt

% Matrixlängen
lengthi = floor(abs(H_max - H_0) / Delta_H + 1);
% lengthj = floor(abs(m_Bat_max - m_Bat_min) / m_Bat_Delta + 1);
% lengthj = length(Abfrage_mot);
% lengthj = length(Abfrage_c_W);
lengthj = length(Abfrage_n_prop);
% lengthj = length(Abfrage_n_bat);
lengthvkg = floor(abs(V_Kg_max - V_Kg_min) / V_Kg_Delta + 1);


%% Initialisierungen

% Initialisierung
l1 = zeros(lengthj,1);
l2 = zeros(lengthj,1);
l3 = zeros(lengthj,1);
l4 = zeros(lengthj,1);
l5 = zeros(lengthj,1);
l6 = zeros(lengthj,1);
l7 = zeros(lengthj,1);
l8 = zeros(lengthj,1);
l9 = zeros(lengthj,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Programmanfang

% j = 1;
% for m_Bat_variabel = m_Bat_min:m_Bat_Delta:m_Bat_max
for j = 1:length(Abfrage_n_prop)
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% allgemeine Parameter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Motor
    motor_name = axi_motor_db{Abfrage_mot(2),1}; % Motorname
    [K_V, I_0, R_i, m_Mot, S_max, I_max] = Motordata('axi_motor_db',motor_name);
    K_V = K_V*2*pi/60;          % Umrechnung in 1/(V*s)
    % R_i = 0.123;            % Innenwiderstand in Ohm
    % K_V = 1400*2*pi/60;     % K_V Wert in 1/(V*s)
    % I_0 = 0.56;             % Leerlaufstrom in Ampere
    % I_max = 30;             % Max Continuous Current
    % m_Mot = 0.0365;         % Motorgewicht in kg
    
    % Propeller
    prop_name = Abfrage_prop{2};    % Propellerbezeichnung
    n_Prop = Abfrage_n_prop(j);             % Anzahl der Propeller
    c_d0 = 0.05;            % Schaetzung des mittleren Nullwiderstandbeiwerts
    a_alpha = 5;            % Anstieg des Auftriebsbeiwerts ueber dem Anstellwinkel (Profil), Schaetzung
    alpha_stall = 10;       % Anstellwinkel, bei dem die Strömung abreisst in Grad, Schaetzung
    
    %Entnahme des Durchmessers und des Pitches aus dem Propellernahmen
    [RPM_map, V_map, T_map, P_map, TAU_map] = Propeller_map(DATA_APC,prop_name);    % Aufbau des Kennfeldes
    D = str2double(prop_name(1:strfind(prop_name,'x')-1));      % Durchmesser extrahieren
    P_75 = prop_name(strfind(prop_name,'x')+1:end);             % Pitch extrahieren
    while isnan(str2double(P_75)) == 1
        P_75(end) = [];
    end
    P_75 = str2double(P_75);                    % Pitch festlegen
    R = D * 0.0254 / 2;                         % Propellerradius in Meter
    F = pi * R^2;                               % Fläche eines Propellers in Quadratmeter
    Theta_75 = atan( 4*P_75 / (3*pi * D) );     % geometrischer Anstellwinkel des Propellers bei 75% des Radius

    % Skalierungsfaktor
    m_ges = (n_Prop * m_Mot)/(4*0.0365/1.06);     % Gesamtmasse
    m_Bat = m_ges * (0.56/1.06);              % Batteriemasse
    m_copter = m_ges * (0.354/1.06);           % Leermasse
    
    % Batterie
    E_Dichte = 890540;      % Energiedichte des LiPos in J/kg
    N_Bat_cell = Abfrage_n_bat(4);         % Anzahl der Batteriezellen in Reihe
    N_Bat_cell_p = 3;       % Anzahl der Batteriezellen parallel
    C_Bat_cell = 3.120;     % Kapazität einer Zelle in Ah
    U_Bat_cell = 4;       % nominale Spannung pro Batteriezelle
    U_Bat_cell_min = 2.85;  % minimale Spannung pro Batteriezelle
    P_Bat_Peukert = 1.00;   % Peukert-Konstante (Schaetzung)
    C_Rate_max = 30;        % maximale C-Rate bezogen auf eine nominale Entladezeit von 1 Stunde
    % m_Bat = 0.56;           % Batteriemasse in kg
    
    % Batteriemassendiskreitisierung
    m_Bat_min = 0.49;
    m_Bat_Delta = 0.01;
    m_Bat_max = 0.55;
    Abfrage_m_Bat = [0.25 0.5 0.75 1 1.5 2];
    
    % Missionsparameter
    m_nutz = 0.0;          % Nutzlast in kg
    
    
    %% Parameter Multicopter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Gesamtsystem
    % m_copter = 0.354;                       % Multicopter Leermasse in kg
    % A_copter = 0.15*0.05 + 0.12*0.02*n_Prop;     % obere Stirnflaeche des Multicopter in m^2
    A_copter = 0.15*0.05 + (D/2*0.0254)*1.2*0.02* n_Prop;     % obere Stirnflaeche des Multicopter in m^2
    A_copter_seitlich = 1.5 * A_copter;     % seitliche Stirnflaeche des Multicopter in m^2
    c_W_copter_oben = Abfrage_c_W(3);                    % Widerstandsbeiwert des Multicopters
    c_W_copter_seitlich = 1 * A_copter_seitlich / A_copter;         % seitlicher Widerstandsbeiwert  des Multicopters
    c_A_copter_max = 0.3;                   % maximaler Auftriebsbeiwert des Multicopters (bei +/-45° Anstellwinkel)
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    U_Bat_nom = N_Bat_cell * U_Bat_cell;        % nominale Batteriespannung
    U_Bat_min = N_Bat_cell * U_Bat_cell_min;    % minimale Batteriespannung   
    C_Bat = E_Dichte * m_Bat / U_Bat_nom;                               % Kapazitaet der Batterie in As
    Cnom = C_Bat/1000;                                            % Nominelle Kapazität
    
    % Umgebung
    H = zeros(lengthi,1);
    rho = zeros(lengthi,1);
    % Multicopter
    alpha = zeros(lengthi,1);
    Theta = zeros(lengthi,1);
    V_Kg = zeros(lengthi,1);
    % Flugzustand_Flaechenflzg = zeros(lengthi,1);
    % V_Flaechenflugzeug = zeros(lengthi,1);
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
        
        % Gesamtlänge der Vektoren: lengthall = lengthgamma + lengthb_vert;
        
        
        
        % Multicopter
        Theta_inter = zeros(lengthvkg, 1);
        alpha_inter = zeros(lengthvkg, 1);
        V_Kg_inter = zeros(lengthvkg,1);
        Bestimmung_V_Kg = zeros(lengthvkg, 1);
        % Propeller
        Omega_inter = zeros(lengthvkg, 1);
        tau_inter = zeros(lengthvkg, 1);
        M_tip_inter = zeros(lengthvkg, 1);
        % Motor
        I_mot_inter = zeros(lengthvkg, 1);
        U_mot_inter = zeros(lengthvkg, 1);
        % ESC
        PWM_inter = zeros(lengthvkg, 1);
        % Batterie
        I_Bat_inter = zeros(lengthvkg, 1);
        U_Bat_inter = zeros(lengthvkg, 1);
        C_Rate_inter = zeros(lengthvkg, 1);
        C_Rest_V_inter = zeros(lengthvkg, 1);
        Delta_C_Bat_inter = zeros(lengthvkg,1);
        i_int_inter = zeros(lengthvkg,1);
        % Gesamtsystem
        Thrust_inter = zeros(lengthvkg, 1);
        eta_ges_inter = zeros(lengthvkg, 1);
        
        z = 1;
        for V_Kg_variabel = V_Kg_min:V_Kg_Delta:V_Kg_max	    % Variation des Bahnneigungswinkel für das Flächenflugzeug
            
            V_Kg_inter(z) = V_Kg_variabel;
            
            % MULTICOPTER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            m = m_copter + m_Bat + m_Mot * n_Prop + m_nutz;                     % Gesamtmasse des Quadrocopters
            t_Flug = Delta_H / V_Kg_inter(z);                                            % Flugzeit
            
            % Aerodynamik
            [Thrust_inter(z),Theta_inter(z),V_A,alpha_inter(z)] = MulticopterAerodynamik(u_Wg,V_Kg_inter(z),gamma_copter,c_W_copter_seitlich,c_W_copter_oben,c_A_copter_max,rho(x),A_copter,m,g);
            
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
                
                [I_Bat_inter(z),U_Bat_inter(z),C_Rate_inter(z),Delta_C_Bat_inter(z),C_Rest_V_inter(z),i_int_inter(z)] = Batterie(Batterie_data,...
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
                    alpha_inter(z) > alpha_stall || M_tip_inter(z) >= 1 || I_Bat_inter(z) <= 0 || eta_ges_inter(z) > 1    % ||  PWM_inter(z) > 1.0
                C_Rest_V_inter(z) = NaN;
                Omega_inter(z) = NaN;
                U_mot_inter(z) = NaN;
                I_mot_inter(z) = NaN;
                I_Bat_inter(z) = NaN;
                PWM_inter(z) = NaN;
                eta_ges_inter(z) = NaN;
                V_Kg_inter(z) = NaN;
                
            end
            
            if x > 1
                if C_Rest_V_inter(z) > C_Rest_V(x-1)
                    C_Rest_V_inter(z) = NaN;
                    Omega_inter(z) = NaN;
                    U_mot_inter(z) = NaN;
                    I_mot_inter(z) = NaN;
                    I_Bat_inter(z) = NaN;
                    PWM_inter(z) = NaN;
                    eta_ges_inter(z) = NaN;
                    V_Kg_inter(z) = NaN;
                end
            end
            
            Bestimmung_V_Kg(z) = Delta_C_Bat_inter(z) * U_Bat_inter(z);        % Berechnung der aufgebrachten Energiemenge
            
            z = z+1;			% Erhöhung der Zählervariablen für die gamma-Schleife
            
        end
        
        
        if isnan(nanmean(Omega_inter)) ~= 1 && isnan(nanmean(I_mot_inter)) ~= 1 &&  isnan(nanmean(U_mot_inter)) ~= 1 && ...     % Wenn alle der Vektoren nicht nur NaN
                isnan(nanmean(PWM_inter)) ~= 1 && isnan(nanmean(C_Rest_V_inter)) ~= 1 && isnan(nanmean(I_Bat_inter)) ~= 1       % enthalten (unfliegbarer Zustand)
            
            % Kriterium für optimalen Steigwinkel
            y = 1;
            while y < length(Bestimmung_V_Kg)                           % Suche und finde optimalen Steiggeschwindigkeit
                A = Bestimmung_V_Kg;
                A = sort(A);                                            % Belegung A mit dem sortierten Werten
                ind_opt = find(Bestimmung_V_Kg == A(y));                % Index mit optimalen Flugzustand
                
                if length(ind_opt) > 1                                  % Wenn die Länge des Index größer als 1 ist ...
                    ind_opt = ind_opt(1);                               % nimm den ersten Eintrag
                end
                
                if ~isnan(Thrust_inter(ind_opt)) && ~isnan(Omega_inter(ind_opt)) && ~isnan(I_mot_inter(ind_opt)) && ~isnan(U_mot_inter(ind_opt)) && ~isnan(PWM_inter(ind_opt)) ...
                        && ~isnan(I_Bat_inter(ind_opt)) && ~isnan(U_Bat_inter(ind_opt)) && ~isnan(C_Rate_inter(ind_opt)) && ~isnan(C_Rest_V_inter(ind_opt))
                    % Falls alle Einträge der Vektoren mit ind_opt real
                    % sind, dann ...
                    break;                                                 % Unterbreche Schleife, falls alle Leistungsparameter physikalisch realistische Werte besitzen
                end
                
                y = y+1;
            end
            
            % Übergabe der Leistungswerte für die optimale Geschwindigkeit und
            % Festlegen für entsprechenden Höhenschritt
            
            % Multicopter
            Theta(x) = Theta_inter(ind_opt);
            alpha(x) = alpha_inter(ind_opt);
            V_Kg(x) = V_Kg_inter(ind_opt);	    % Bestimmung opt. Stgeschw. und speichern in Vektor für jeden Höhenabschnitt
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
            % Multicopter
            Theta(x) = NaN;
            alpha(x) = NaN;
            V_Kg(x) = NaN;
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
        
        
        H(x) = H_oben;			% Speichern der Höhe im Vektor
        x = x+1;				% Erhöhung der Zählervariablen für die Höhen-Schleife
        disp([num2str((x-1)*100/lengthi) ' %']);
        
    end
    
    % Darstellung der Ergenisse in Diagrammen
    figure(figure_ges)
    
    %     subplot(521), l1(j) = plot(H,C_Rest_V*100,'LineWidth',2); l1_Info{j} = ['m_{Bat} = ' num2str(m_Bat) ' kg']; grid on, hold on
    %     subplot(522), l2(j) = plot(H,Omega/(2*pi)*60,'LineWidth',2); l2_Info{j} = ['m_{Bat} = ' num2str(m_Bat) ' kg']; grid on, hold on
    %     subplot(523), l3(j) = plot(H,I_mot,'LineWidth',2); l3_Info{j} = ['m_{Bat} = ' num2str(m_Bat) ' kg']; grid on, hold on
    %     subplot(524), l4(j) = plot(H,U_mot,'LineWidth',2); l4_Info{j} = ['m_{Bat} = ' num2str(m_Bat) ' kg']; grid on, hold on
    %     subplot(525), l5(j) = plot(H,I_Bat,'LineWidth',2); l5_Info{j} = ['m_{Bat} = ' num2str(m_Bat) ' kg']; grid on, hold on
    %     H2 = [0;H];
    %     subplot(526), l6(j) = plot(H2,U_Bat,'LineWidth',2); l6_Info{j} = ['m_{Bat} = ' num2str(m_Bat) ' kg']; grid on, hold on
    %     subplot(527), l7(j) = plot(H,PWM*100,'LineWidth',2); l7_Info{j} = ['m_{Bat} = ' num2str(m_Bat) ' kg']; grid on, hold on
    %     subplot(528), l8(j) = plot(H,eta_ges*100,'LineWidth',2); l8_Info{j} = ['m_{Bat} = ' num2str(m_Bat) ' kg']; grid on, hold on
    %     subplot(529), l9(j) = plot(H,V_Kg,'LineWidth',2); l9_Info{j} = ['m_{Bat} = ' num2str(m_Bat) ' kg']; grid on, hold on
    
    subplot(521), l1(j) = stairs(H,C_Rest_V*100,'LineWidth',1); grid on, hold on
    subplot(522), l2(j) = stairs(H,Omega/(2*pi)*60,'LineWidth',1); grid on, hold on
    subplot(523), l3(j) = stairs(H,I_mot,'LineWidth',1); grid on, hold on
    subplot(524), l4(j) = stairs(H,U_mot,'LineWidth',1); grid on, hold on
    subplot(525), l5(j) = stairs(H,I_Bat,'LineWidth',1); grid on, hold on
    H2 = [0;H];
    subplot(526), l6(j) = stairs(H2,U_Bat,'LineWidth',1); grid on, hold on
    subplot(527), l7(j) = stairs(H,PWM*100,'LineWidth',1); grid on, hold on
    subplot(528), l8(j) = stairs(H,eta_ges*100,'LineWidth',1); grid on, hold on
    subplot(529), l9(j) = stairs(H,V_Kg,'LineWidth',1); 
    l9_Info{j} = ['n_{Prop} =' num2str(Abfrage_n_prop(j))]; grid on, hold on
    
    
%     'm_{Mot}= ' num2str(m_Mot) ', K_V= ' num2str(K_V*60/(2*pi)) ', Prop= ' prop_name ', n_{Bat,cell}= ' num2str(N_Bat_cell)
%   ['n_{Prop} =' num2str(Abfrage_n_prop(j))]  
    
    
    
%     j = j + 1;
    %% Spielereien
    disp([num2str((j)*100/lengthj) ' %']); %(Abfrage_m_Bat)) ' %']);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot

% Darstellung der Ergenisse in Diagrammen
figure(figure_ges)

subplot(521), title('Restladung'), xlabel('Höhe [m]'),
ylabel('C_{Bat,Rest} [%]'),% legend([l1], l1_Info,'Location','northeastoutside')
subplot(522), title('Drehzahl'), xlabel('Höhe [m]'),
ylabel('\Omega [RPM]'), %legend([l2], l2_Info)
subplot(523), title('Motorstrom'), xlabel('Höhe [m]'),
ylabel('I_{Mot} [A]'),% legend([l3], l3_Info)
subplot(524), title('Motorspannung'), xlabel('Höhe [m]'),
ylabel('U_{mot} [V]'),% legend([l4], l4_Info)
subplot(525), title('Batteriestrom'), xlabel('Höhe [m]'),
ylabel('I_{Bat} [A]'), %legend([l5], l5_Info)
H2 = [0;H];
subplot(526), title('Batteriespannung'), xlabel('Höhe [m]'),
ylabel('U_{Bat} [A]'), %legend([l6], l6_Info)
subplot(527), title('Pulsweitenmodulation'), xlabel('Höhe [m]'),
ylabel('PWM [%]'),% legend([l7], l7_Info)
subplot(528), title('Gesamtwirkungsgrad'), xlabel('Höhe [m]'),
ylabel('\eta_{ges} [%]'), %legend([l8], l8_Info)
subplot(529), title('Bahngeschwindigkeit'), xlabel('Höhe [m]'),
ylabel('V_{Kg} [m/s]'),
lgd = legend([l9], l9_Info, 'Location','bestoutside'); title(lgd,'Größenskalierung')  
% 
% newPosition = [10 20 0.2 0.2];
% newUnits = 'centimeters';
% set(hL,'Position', newPosition,'Units', newUnits);

% Anpassung und Abspeichern der Diagramme
PaperSizeX = 21;
PaperSizeY = 29.7;

fig = gcf;
set(gcf,'PaperUnits','centimeters', 'PaperPosition', [-2 -2.6 24.65 34.45]);%[-1.75 -2.6 24.65 34.45]);
set(gcf,'Units','centimeters', 'PaperSize', [PaperSizeX PaperSizeY]);
saveas(gcf,Dateiname, 'pdf');