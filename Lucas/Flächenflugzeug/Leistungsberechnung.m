%% Leistungsberechnung f�r Flugsysteme %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% *************************************************************************


%% Initialisierungen %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Batterie

U_Bat_nom = N_Bat_cell * U_Bat_cell;        % nominale Batteriespannung
U_Bat_min = N_Bat_cell * U_Bat_cell_min;    % minimale Batteriespannung
C_Bat = E_Dichte * m_Bat / U_Bat_nom;       % Kapazitaet der Batterie in As
% C_Bat = N_Bat_cell_p*C_Bat_cell*3600;       % Kapazit�t in As
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



%% Umgebungsparameter

T_11 = T_0 - 0.0065 *(11000-H_0);                   % T in 11000m H�he
rho_11 = rho_0 * (1 - 0.0065*(11000/T_0))^4.256;    % Dichte in 11000m H�he
p_11 = p_0 * (1 - 0.0065*(11000/T_0))^5.256;        % Druck in 11000m H�he



%% Intialisierung der Matrizen f�r jeden H�henabschnitt

lengthi = floor(abs(H_max - H_0) / Delta_H + 1);
lengthgamma = floor(abs(gamma_max - gamma_min) / gamma_Delta + 1);
v_ver_max = ceil(V_stern);
lengthvert = floor(abs(v_ver_max - v_vert_min) / v_vert_Delta + 1);
% lengthj = length(Abfrage_E);
% lengthj = length(Abfrage_V);
lengthj = length(Abfrage_f_P);


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
l10 = zeros(lengthj,1);

% j = 1;
for j = 1:length(Abfrage_f_P)
    
%     E_stern = Abfrage_E(j);
    E = E_stern;
%     V_stern = Abfrage_V(j)/3.6;
    
    m_Mot_Quad = 0.0365;
    n_Prop_Quad = 4;
    
%     m = m_copter + n_Prop_Quad*m_Mot_Quad + m_nutz + m_Bat;
    
    f_p = Abfrage_f_P(j);                                    % Penalty-Faktor f�r das Strukturgewicht des Flugzeugs
    
    m_flugzeug = f_p * m_copter;
    m_Bat = 0.56;
    m_Bat = m_Bat+((n_Prop_Quad*m_Mot_Quad)-m_Mot*n_Prop) + (1-f_p) * m_copter;


    [RPM_map, V_map, T_map, P_map, TAU_map] = Propeller_map(DATA_APC,prop_name);    % Aufbau des Kennfeldes
    
    % Umgebung
    H = zeros(lengthi,1);
    rho = zeros(lengthi,1);
    % Multicopter
    alpha = zeros(lengthi,1);
    Theta = zeros(lengthi,1);
    % Fl�chenflugzeug
    gamma_Flaechenflzg = zeros(lengthi,1);
    Flugzustand_Flaechenflzg = zeros(lengthi,1);
    V_Flaechenflugzeug = zeros(lengthi,1);
    % Propeller
    Thrust = zeros(lengthi,1);
    Omega = zeros(lengthi,1);
    tau = zeros(lengthi,1);
    M_tip = zeros(lengthi,1);
    eta_prop = zeros(lengthi,1);
    % Motor
    U_mot = zeros(lengthi,1);
    I_mot = zeros(lengthi,1);
    eta_mot = zeros(lengthi,1);
    % ESC
    PWM = zeros(lengthi,1);
    eta_PWM = zeros(lengthi,1);
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
    eta_ges = zeros(lengthi,1);
    
    
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
        
        % Gesamtl�nge der Vektoren: lengthall = lengthgamma + lengthb_vert;
        
        
        % Multicopter
        Theta_inter = zeros(lengthgamma+lengthvert, 1);
        alpha_inter = zeros(lengthgamma+lengthvert, 1);
        % Flaechenflugzeug
        Flugzustand_Flaechenflzg_inter = zeros(lengthgamma+lengthvert,1);
        gamma_Flaechenflzg_inter = zeros(lengthgamma+lengthvert,1);
        Bestimmung_gamma = zeros(lengthgamma, 1);
        V_Flaechenflugzeug_inter = zeros(lengthgamma+lengthvert, 1);
        % Propeller
        Omega_inter = zeros(lengthgamma+lengthvert, 1);
        tau_inter = zeros(lengthgamma+lengthvert, 1);
        M_tip_inter = zeros(lengthgamma+lengthvert, 1);
        eta_prop_inter = zeros(lengthgamma+lengthvert,1);
        % Motor
        I_mot_inter = zeros(lengthgamma+lengthvert, 1);
        U_mot_inter = zeros(lengthgamma+lengthvert, 1);
        eta_mot_inter = zeros(lengthgamma+lengthvert,1);
        % ESC
        PWM_inter = zeros(lengthgamma+lengthvert, 1);
        eta_PWM_inter = zeros(lengthgamma+lengthvert,1);
        % Batterie
        I_Bat_inter = zeros(lengthgamma+lengthvert, 1);
        U_Bat_inter = zeros(lengthgamma+lengthvert, 1);
        C_Rate_inter = zeros(lengthgamma+lengthvert, 1);
        C_Rest_V_inter = zeros(lengthgamma+lengthvert, 1);
        Delta_C_Bat_inter = zeros(lengthgamma+lengthvert,1);
        i_int_inter = zeros(lengthgamma+lengthvert,1);
        % Gesamtsystem
        Thrust_inter = zeros(lengthgamma+lengthvert, 1);
        eta_ges_inter = zeros(lengthgamma+lengthvert, 1);
        
        b = 0;                          % Kontrollfaktor
        V_vert = 0;                     % Vertikalgeschwindigkeit Z�hlfaktor
        
        z = 1;
        
        for gamma_variabel = gamma_min:gamma_Delta:(gamma_max+v_ver_max)	    % Variation des Bahnneigungswinkel f�r das Fl�chenflugzeug
            
            gamma_Flaechenflzg_inter(z) = gamma_variabel;
            
            if Abfrage_Flugsystem == 1
                
                % MULTICOPTER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                m = m_copter + m_Bat + m_Mot * n_Prop + m_nutz;                     % Gesamtmasse des Quadrocopters
                t_Flug = Delta_H / V_Kg;                                            % Flugzeit
                
                % Aerodynamik
                [Thrust_inter(z),Theta_inter(z),V_A,alpha_inter(z)] = MulticopterAerodynamik(u_Wg,V_Kg,gamma_copter,c_W_copter_seitlich,c_W_copter_oben,c_A_copter_max,rho(x),A_copter,m,g);
                
                
            else
                
                
                % FL�CHENFLUGZEUG %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                m = m_flugzeug + m_Bat + m_Mot * n_Prop + m_nutz;                   % Gesamtmasse des Fl�chenflugzeugs
                
                % Hier fehlt noch Behandlung des Zustandes, wenn kein
                % Flugzustand == 2 erreicht wird
                
                
                % Aerodynamik
                if b == 0                      % Solange Flugzustand = 2 noch nicht erreicht ist, Berechne die Aerodynamik nach ...
                    [Thrust_inter(z),V_A,Flugzustand_Flaechenflzg_inter(z)] = FlaechenflugzeugAerodynamik(m,g,E_stern,V_stern,rho_stern,E,gamma_variabel,rho(x));
                    V_H = V_A * sind(gamma_variabel);         % Bestimmung der Vertikalgeschwindigkeit
                    t_Flug = Delta_H / V_H;                   % Berechne aus Vertikalgeschwindigkeit und H�he die Flugzeit
                    V_Flaechenflugzeug_inter(z) = V_A;        % Speicher die absolute Fluggeschwindigkeit
                end
                
                if Flugzustand_Flaechenflzg_inter(z) == 2 || b == 1 || gamma_variabel >= 90   % Wenn Flugzustand = 2 erreicht wurde
                    % oder der Steigwinkel bereits 90� ist, berechne die Aerodynamik nach
                    b = 1;                                     % Setze den Kontrollfaktor gleich 1
                    V_vert = V_vert + v_vert_Delta;            % Iteriere die Vertikalgeschwindigkeit
                    V_A = V_vert;                              % Vertikalgeschw. = Fluggeschw.
                    gamma_Flaechenflzg_inter(z) = 90;          % Der Bahnneigungswinkel ist 90�
                    [Thrust_inter(z)] = Vertikalflug(m,g,E_stern,V_stern,rho_stern,V_A,rho(x));
                    t_Flug = Delta_H / V_A;                    % Berechne aus Vertikalgeschwindigkeit und H�he die Flugzeit
                    V_Flaechenflugzeug_inter(z) = V_A;         % Speicher die absolute Fluggeschwindigkeit
                end
                
                
                alpha_inter(z) = - 90 * pi/180;                % Die Rotorebene f�r ein Fl�chenflugzeug ist immer 90�
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
                
                
                % Wie gro� ist die Blattspitzengeschwindigkeit?
                M_tip_inter(z) = (Omega_inter(z) * R)/a;                       % Blattspitzengeschwindigkeit in Ma
                
                
                % Motorzustand berechnen
                [U_mot_inter(z),I_mot_inter(z),eta_mot_inter(z)] = Motor(tau_inter(z),K_V,I_0,R_i,Omega_inter(z));
                
                
                % Zustand der Motorregler berechnen
                [PWM_inter(z),eta_PWM_inter(z)] = ESC(U_mot_inter(z),U_Bat(x));         % <-- hier U_bat_inter
                
                
                % Batteriezustand berechnen
                Delta_C_Bat_inter(z) = Delta_C_Bat(x);                         % Anpassung der Batteriekapazit�t
                i_int_inter(z) = i_int(x);                                     % �bergabe des Integrals der Spannung vom letzten Schritt
                %            [I_Bat_inter(z),C_Rate_inter(z),Delta_C_Bat_inter(z),C_Rest_V_inter(z)] = Batterie(PWM_inter(z),eta_PWM,I_mot_inter(z),n_Prop,C_Bat,P_Bat_Peukert,Delta_C_Bat_inter(z),t_Flug);
                
                [I_Bat_inter(z),U_Bat_inter(z),C_Rate_inter(z),Delta_C_Bat_inter(z),C_Rest_V_inter(z),i_int_inter(z)] = Batterie_neu(Batterie_data,...
                    Cnom,PWM_inter(z),eta_PWM_inter(z),n_Prop,i_int_inter(z),U_Bat_inter(z),C_Bat,Delta_C_Bat_inter(z),I_mot_inter(z),N_Bat_cell,P_Bat_Peukert,t_Flug);
                
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
                eta_prop_inter(z) = (Thrust_inter(z) * (mu_z + vi))/(tau_inter(z) * Omega_inter(z));
                
                eta_ges_inter(z) = (n_Prop * Thrust_inter(z) * (mu_z + vi))/(I_Bat_inter(z) * U_Bat_inter(z));         % Leistung, die in Schub umgesetzt wird im Verh�ltnis zur aufgebrachten Leistung
                
                
            end
            
            % Wenn Grenzen ueberschritten werden, Resultate entfernen
            alpha_inter(z) = alpha_inter(z)*180/pi;
            
            if C_Rest_V_inter(z) < 0.0 || U_mot_inter(z) > U_Bat(x) || U_mot_inter(z) <= 0 || C_Rate_inter(z) > C_Rate_max || I_mot_inter(z) > I_max || ...
                    alpha_inter(z) > alpha_stall || M_tip_inter(z) >= 1 || I_Bat_inter(z) <= 0  || eta_ges_inter(z) > 1  %|| C_Rest_V_inter(z) > C_Rest_V(x-1) ||  PWM_inter(z) > 1.0
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
            
            if x > 1     % Entferne alle Ergebnisse, wobei die Restladung im Vergleich zum vorheigen Schritt steigt
                if C_Rest_V_inter(z) > C_Rest_V(x-1)
                    C_Rest_V_inter(z) = NaN;
                    Omega_inter(z) = NaN;
                    U_mot_inter(z) = NaN;
                    I_mot_inter(z) = NaN;
                    I_Bat_inter(z) = NaN;
                    PWM_inter(z) = NaN;
                    eta_ges_inter(z) = NaN;
                    V_Flaechenflugzeug_inter(z) = NaN;
                end
            end
            
            % Muss der Steigwinkel variiert werden?
            
            if Abfrage_Flugsystem == 1
                break;								% Abbruch f�r den Fall eines Multicopters um for-Schleife zu verlassen, da gamma vorgegeben
            else
                Bestimmung_gamma(z) = Delta_C_Bat_inter(z) * U_Bat_inter(z);        % Berechnung der aufgebrachten Energiemenge
            end
            
            
            z = z+1;			% Erh�hung der Z�hlervariablen f�r die gamma-Schleife
            
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
            eta_prop(x) = eta_prop_inter(1);
            % Motor
            I_mot(x) = I_mot_inter(1);
            U_mot(x) = U_mot_inter(1);
            eta_mot(x) = eta_mot_inter(1);
            % ESC
            PWM(x) = PWM_inter(1);
            eta_PWM(x) = eta_PWM_inter(1);
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
                
                % Kriterium f�r optimalen Steigwinkel
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
                
                % �bergabe der Leistungswerte f�r die optimale Geschwindigkeit und
                % Festlegen f�r entsprechenden H�henschritt
                
                % Fl�chenflugzeug
                Theta(x) = Theta_inter(ind_opt);
                alpha(x) = alpha_inter(ind_opt);
                gamma_Flaechenflzg(x) = gamma_Flaechenflzg_inter(ind_opt);	    % Bestimmung opt. Stgeschw. und speichern in Vektor f�r jeden H�henabschnitt
                Flugzustand_Flaechenflzg(x) = Flugzustand_Flaechenflzg_inter(ind_opt);
                V_Flaechenflugzeug(x) = V_Flaechenflugzeug_inter(ind_opt);
                % Propeller
                Thrust(x) = Thrust_inter(ind_opt);
                Omega(x) = Omega_inter(ind_opt);
                tau(x) = tau_inter(ind_opt);
                M_tip(x) = M_tip_inter(ind_opt);
                eta_prop(x) = eta_prop_inter(ind_opt);
                % Motor
                I_mot(x) = I_mot_inter(ind_opt);
                U_mot(x) = U_mot_inter(ind_opt);
                eta_mot(x) = eta_mot_inter(ind_opt);
                % ESC
                PWM(x) = PWM_inter(ind_opt);
                eta_PWM(x) = eta_PWM_inter(ind_opt);
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
                % Fl�chenflugzeug
                Theta(x) = NaN;
                alpha(x) = NaN;
                gamma_Flaechenflzg(x) = NaN;
                V_Flaechenflugzeug(x) = NaN;
                % Propeller
                Thrust(x) = NaN;
                Omega(x) = NaN;
                tau(x) = NaN;
                M_tip(x) = NaN;
                eta_prop(x) = NaN;
                % Motor
                I_mot(x) = NaN;
                U_mot(x) = NaN;
                eta_mot(x) = NaN;
                % ESC
                PWM(x) = NaN;
                eta_PWM(x) = NaN;
                % Batterie
                I_Bat(x) = NaN;
                U_Bat(x+1) = NaN;
                C_Rate(x) = NaN;
                C_Rest_V(x) = NaN;
                % Gesamtsystem
                eta_ges(x) = NaN;
                
            end
            
        end
        
        
        H(x) = H_oben;			% Speichern der H�he im Vektor
        x = x+1;				% Erh�hung der Z�hlervariablen f�r die H�hen-Schleife
        
        disp([num2str((x-1)*100/lengthi) ' %']);
    end
    figure(figure_ges)
   
    subplot(521), l1(j) = stairs(H,C_Rest_V*100,'LineWidth',1); grid on, hold on
    subplot(522), l2(j) = stairs(H,Omega/(2*pi)*60,'LineWidth',1); grid on, hold on
    subplot(523), l3(j) = stairs(H,I_mot,'LineWidth',1); grid on, hold on
    subplot(524), l4(j) = stairs(H,U_mot,'LineWidth',1); grid on, hold on
    subplot(525), l5(j) = stairs(H,I_Bat,'LineWidth',1); grid on, hold on
    H2 = [0;H];
    subplot(526), l6(j) = stairs(H2,U_Bat,'LineWidth',1); grid on, hold on
    subplot(527), l7(j) = stairs(H,PWM*100,'LineWidth',1); grid on, hold on
    subplot(528), l8(j) = stairs(H,eta_ges*100,'LineWidth',1); grid on, hold on
    subplot(529), l9(j) = stairs(H,gamma_Flaechenflzg,'LineWidth',1); grid on, hold on
    subplot(5,2,10), l10(j) = stairs(H,V_Flaechenflugzeug,'LineWidth',1); l10_Info{j} = ['f_P = ' num2str(Abfrage_f_P(j))]; grid on, hold on
    
    %% Spielereien
    disp([num2str((j)*100/lengthj) ' %']); %(Abfrage_m_Bat)) ' %']);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot

% Darstellung der Ergenisse in Diagrammen  
figure(figure_ges)

subplot(521), title('Restladung'), xlabel('H�he [m]'),
    ylabel('C_{Bat,Rest} [%]'),% legend([l1], l1_Info,'Location','northeastoutside')
subplot(522), title('Drehzahl'), xlabel('H�he [m]'),
    ylabel('Drehzahl [RPM]'), %legend([l2], l2_Info)
subplot(523), title('Motorstrom'), xlabel('H�he [m]'),
    ylabel('I_{Mot} [A]'),% legend([l3], l3_Info)
subplot(524), title('Motorspannung'), xlabel('H�he [m]'),
    ylabel('U_{mot} [V]'),% legend([l4], l4_Info)
subplot(525), title('Batteriestrom'), xlabel('H�he [m]'),
    ylabel('I_{Bat} [A]'), %legend([l5], l5_Info)
H2 = [0;H];
subplot(526), title('Batteriespannung'), xlabel('H�he [m]'),
    ylabel('U_{Bat} [A]'), %legend([l6], l6_Info)
subplot(527), title('Pulsweitenmodulation'), xlabel('H�he [m]'),
    ylabel('PWM [%]'),% legend([l7], l7_Info)
subplot(528), title('Gesamtwirkungsgrad'), xlabel('H�he [m]'),
    ylabel('eta_{ges} [%]'), %legend([l8], l8_Info)
subplot(529), title('Bahnneigungswinkel'), xlabel('H�he [m]'),
    ylabel('gamma [�]'), 
subplot(5,2,10), title('Bahngeschwindigkeit'), xlabel('H�he [m]'),
    ylabel('V_{Kg} [m/s]'),    
lgd = legend([l10], l10_Info, 'Location','bestoutside'); title(lgd,'Penalty-Faktor')
    

ImageSizeX = 21;
ImageSizeY = 29.7;
% figure(figure_ges)
% set(gcf,'PaperUnits','centimeters', 'PaperPosition', [0 0 ImageSizeX ImageSizeY]);
% set(gcf,'Units','centimeters', 'PaperSize', [ImageSizeX ImageSizeY]);
% saveas(gcf,Dateiname, 'pdf');

fig = gcf;
fig.PaperPositionMode = 'auto';
set(fig,'PaperUnits','centimeters', 'PaperPosition', [0 0 ImageSizeX ImageSizeY]);
fig_pos = fig.PaperPosition;
fig.PaperSize = [ImageSizeX ImageSizeY];
print(fig,Dateiname,'-dpdf')


