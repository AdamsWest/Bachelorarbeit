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

R = D * 0.0254 / 2;                         % Propellerradius in Meter
F = pi * R^2;                               % Fläche eines Propellers in Quadratmeter

DATA = evalin('base','DATA_APC');

[l,w] = size(DATA_APC);
clear w;
% Aufsuchen und finden aller Propellernamen, die den gleichen Durchmesser
% wie der gesuchte haben. Dabei werden alle anderen Propeller entfernt

for w = l:-1:1
    a = DATA{w,1};                                  % Speichern des Propellernames unter der Variablen a
    diameter = str2double(a(1:strfind(a,'x')-1));   % Extrahieren des Durchmessers aus dem Propellernamen
    if diameter ~= D                                % Entfernen der Zeile mit einem anderen Propellerdurchmesser
        DATA(w,:) = [];
    end
end

% Entnahme jedes einzelnen Propellers
counter = 0;
ind = 1;
len = length(DATA);
d = 1;
while ind < len
    
    prop_name = DATA{ind,1};                                % Übergabe des Propellernamen
    pitch_all = prop_name(strfind(prop_name,'x')+1:end);    % Übergabe des Pitch
    ind = find(strcmp(DATA(:,1),prop_name));                % Finde und speicher alle Indizes mit diesem Propellernamen
    
    ind_unten = min(ind);                                   % Unteres Ende des Indexvektors
    ind_oben = max(ind);                                    % Oberes ----------"-----------
    
    
    % nicht passende Propeller werden entfernt (-E, 3- and 4-blade, -SF, etc.)
    if isnan(str2double(pitch_all)) == 1                    % Wenn der Propeller eine Zusatzbezeichnung hat ...
        
        for d = ind_oben:-1:ind_unten                       % Gehe die Datenbank durch 
            DATA(d,:) = [];                                 % und entferne diesen Propeller
            len = len - 1;                                  % Für jede gelöschte Zeile verkürze die Länge um 1
        end
        
        ind = min(ind);                                     % Übergabe des neuen Startpunkts für den nächsten Propeller
        
    else
        
        ind = max(ind) + 1;                                 % Übergabe des neuen Startpunkts für den nächsten Propeller
        
        [RPM, V, T, P, Tau] = Propeller_map(DATA,prop_name);% Transformiere das Kennfeld
        counter = counter + 1;                              % Erhöhe die Anzahl der verwendbaren Propeller um 1
        
        % Die Daten für RPM, V, T, usw. werden mit fortlaufenden Nummern Vektoren
        % zugewiesen
        assignin ('base',['prop_name' num2str(counter)], prop_name);
        assignin ('base',['pitch_' num2str(counter)], pitch_all);
        assignin ('base',['RPM_' num2str(counter)], RPM);
        assignin ('base',['V_' num2str(counter)], V);
        assignin ('base',['T_' num2str(counter)], T);
        assignin ('base',['P_' num2str(counter)], P);
        assignin ('base',['Tau_' num2str(counter)], Tau);
    end
end



%% Umgebungsparameter

Temp_11 = T_0 - 0.0065 *(11000-H_0);                   % T in 11000m Höhe
rho_11 = rho_0 * (1 - 0.0065*(11000/T_0))^4.256;    % Dichte in 11000m Höhe
p_11 = p_0 * (1 - 0.0065*(11000/T_0))^5.256;        % Druck in 11000m Höhe



%% Intialisierung der Matrizen für jeden Höhenabschnitt

% Matrixlängen
lengthi = floor(abs(H_max - H_0) / Delta_H + 1);
lengthvkg = floor(abs(V_Kg_max - V_Kg_min) / V_Kg_Delta + 1);

% Umgebung
H = zeros(lengthi,1);
rho = zeros(lengthi,1);
% Multicopter
alpha = zeros(lengthi,1);
Theta = zeros(lengthi,1);
V_Kg_vpp = zeros(lengthi,1);
t_Flug_vpp = zeros(lengthi+1,1);
% Propeller
Thrust = zeros(lengthi,1);
Omega_vpp = zeros(lengthi,1);
tau = zeros(lengthi,1);
M_tip = zeros(lengthi,1);
pitch_vpp = zeros(lengthi,1);
% Motor
U_mot_vpp = zeros(lengthi,1);
I_mot_vpp = zeros(lengthi,1);
% ESC
PWM_vpp = zeros(lengthi,1);
% Batterie
C_Rate = zeros(lengthi,1);
C_Rest_V_vpp = zeros(lengthi,1);
I_Bat_vpp = zeros(lengthi,1);
U_Bat_vpp = zeros(lengthi+1,1);
U_Bat_vpp(1) = U_Bat_nom;                               % Startspannung ist die Nominalspannung
P_Bat = zeros(lengthi,1);
Delta_C_Bat = zeros(lengthi+1,1);
i_int = zeros(lengthi+1,1);
% Gesamtsystem
eta_prop = zeros(lengthi,1);
eta_ges_vpp = zeros(lengthi,1);


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
        rho_oben = rho_11 * exp(-g/(287*Temp_11)*(H_oben-11000));
        
        p_unten = p_0 *(1-0.0065*(H_unten/T_0))^5.256;
        p_oben = p_11 * exp(-g/(287.1*Temp_11)*(H_oben-11000));
    else
        rho_unten = rho_11 * exp(-g/(287*Temp_11)*(H_unten-11000));
        rho_oben = rho_11 * exp(-g/(287*Temp_11)*(H_oben-11000));
        
        p_unten = p_11 * exp(-g/(287.1*Temp_11)*(H_unten-11000));
        p_oben = p_11 * exp(-g/(287.1*Temp_11)*(H_oben-11000));
    end
    
    if H_mitte <= 11000                                             % mittlere Temperatur im Intervall
        T = T_0 - 0.0065 * H_mitte;
    else
        T = Temp_11;
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
    
    
    %     T_map = T_map * rho(x)/rho_1;                                   % Anpassung des Schubkennfeldes an die sich ändernde Dichte
    %     P_map = P_map * rho(x)/rho_1;                                   % Anpassung des Leistungskennfeldes an die sich ändernde Dichte
    %     TAU_map = TAU_map * rho(x)/rho_1;                               % Anpassung des Drehmomentkennfeldes an die sich ändernde Dichte
    
    
    %% Beginn der Leistungsberechnung für Flugsysteme
    
    % Initialisierung des Auswahlvektors für den Steigflug
    
    % Gesamtlänge der Vektoren: lengthall = lengthgamma + lengthb_vert;
    
    
    % Multicopter
    Theta_inter = zeros(lengthvkg, 1);
    alpha_inter = zeros(lengthvkg, 1);
    V_Kg_inter = zeros(lengthvkg,1);
    Bestimmung_V_Kg = zeros(lengthvkg, 1);
    t_Flug_inter = zeros(lengthvkg,1);
    % Propeller
    Omega_inter = zeros(lengthvkg, 1);
    tau_inter = zeros(lengthvkg, 1);
    M_tip_inter = zeros(lengthvkg, 1);
    pitch_inter = zeros(lengthvkg, 1);
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
        t_Flug_inter(z) = Delta_H / V_Kg_inter(z);                                            % Flugzeit
        
        % Initialisierungen für den Verstellpropeller
        % Multicopter
        Theta_vprop = zeros(counter, 1);
        alpha_vprop = zeros(counter, 1);
        V_Kg_vprop = zeros(counter, 1);
        Bestimmung_vprop = zeros(counter, 1);
        % Propeller
        Omega_vprop = zeros(counter, 1);
        tau_vprop = zeros(counter, 1);
        M_tip_vprop = zeros(counter, 1);
        pitch_vprop = zeros(counter,1);
        % Motor
        I_mot_vprop = zeros(counter, 1);
        U_mot_vprop = zeros(counter, 1);
        % ESC
        PWM_vprop = zeros(counter, 1);
        % Batterie
        I_Bat_vprop = zeros(counter, 1);
        U_Bat_vprop = zeros(counter, 1);
        C_Rate_vprop = zeros(counter, 1);
        C_Rest_V_vprop = zeros(counter, 1);
        Delta_C_Bat_vprop = zeros(counter, 1);
        i_int_vprop = zeros(counter, 1);
        % Gesamtsystem
        Thrust_vprop = zeros(counter, 1);
        eta_ges_vprop = zeros(counter, 1);
       
        
        for n = 1:counter
            
            RPM_map = evalin('base',['RPM_' num2str(n)]);
            V_map = evalin('base',['V_' num2str(n)]);
            T_map = evalin('base',['T_' num2str(n)]) * rho(x)/rho_0;
            P_map = evalin('base',['P_' num2str(n)]) * rho(x)/rho_0;
            TAU_map = evalin('base',['Tau_' num2str(n)]) * rho(x)/rho_0;          
            
            V_Kg_vprop(n) = V_Kg_inter(z);
            
            % Aerodynamik
            [Thrust_vprop(n),Theta_vprop(n),V_A,alpha_vprop(n)] = MulticopterAerodynamik(u_Wg,V_Kg_vprop(n),gamma_copter,c_W_copter_seitlich,c_W_copter_oben,c_A_copter_max,rho(x),A_copter,m,g);
                      
            Thrust_vprop(n) = Thrust_vprop(n) / n_Prop;                        % Schub auf n Propeller verteilen
            
            if Thrust_vprop(n) > max(max(T_map))                               % wenn Schub zu gross (Ergebnis verwerfen)
                
                Omega_vprop(n) = NaN;
                I_mot_vprop(n) = NaN;
                C_Rate_vprop(n) = NaN;
                C_Rest_V_vprop(n) = NaN;
                
                
            else

                % Drehzahl und Drehmoment bestimmen
                [Omega_vprop(n),tau_vprop(n)] = Propeller(V_A, alpha_vprop(n), Thrust_vprop(n), RPM_map, V_map, T_map, TAU_map);
   
                % Pitch
                zw = str2double(evalin('base',['pitch_' num2str(n)]));
                pitch_vprop(n) = zw;
                
                % Wie groß ist die Blattspitzengeschwindigkeit?
                M_tip_vprop(n) = (Omega_vprop(n) * R)/a;                       % Blattspitzengeschwindigkeit in Ma
                
                
                % Motorzustand berechnen
                [U_mot_vprop(n),I_mot_vprop(n)] = Motor(tau_vprop(n),K_V,I_0,R_i,Omega_vprop(n));
                
                
                % Zustand der Motorregler berechnen
                [PWM_vprop(n),eta_PWM] = ESC(U_mot_vprop(n),U_Bat_vpp(x));         % <-- hier U_bat_inter

                
                % Batteriezustand berechnen
                Delta_C_Bat_vprop(n) = Delta_C_Bat(x);                         % Anpassung der Batteriekapazität
                i_int_vprop(n) = i_int(x);                                     % Übergabe des Integrals der Spannung vom letzten Schritt
                
                [I_Bat_vprop(n),U_Bat_vprop(n),C_Rate_vprop(n),Delta_C_Bat_vprop(n),C_Rest_V_vprop(n),i_int_vprop(n)] = Batterie(Batterie_data,...
                    Cnom,PWM_vprop(n),eta_PWM,n_Prop,i_int_vprop(n),U_Bat_vprop(n),C_Bat,Delta_C_Bat_vprop(n),I_mot_vprop(n),N_Bat_cell,P_Bat_Peukert,t_Flug_inter(z));
                
                
                %% Gesamtwirkungsgrad
                
                % Berechnung der induzierten Geschwindigkeiten nach van der Wall
                % (Grundlagen der Hubschrauber-Aerodynamik) (2015) S. 153
                
                vi0 = sqrt(m*g / ( 2*rho(x)*F*n_Prop ) );                          % induzierte Geschwindigkeit im Schwebeflug v_i0
                v = vi0;
                mu_z = -V_A*sin(alpha_vprop(n));                                          % Geschwindigkeit durch die Rotorebene
                mu = V_A*cos(alpha_vprop(n));                                             % Geschwindigkeit entlang Rotorebene
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
                
                eta_ges_vprop(n) = (n_Prop * Thrust_vprop(n) * (mu_z + vi))/(I_Bat_vprop(n) * U_Bat_vprop(n));         % Leistung, die in Schub umgesetzt wird im Verhältnis zur aufgebrachten Leistung

                
            end
            
            %% Wenn Grenzen ueberschritten werden, Resultate entfernen
            alpha_vprop(n) = alpha_vprop(n)*180/pi;
            
            if C_Rest_V_vprop(n) < 0.0 || U_mot_vprop(n) > U_Bat_vpp(x) || U_mot_vprop(n) <= 0 || C_Rate_vprop(n) > C_Rate_max || I_mot_vprop(n) > I_max || ...
                    alpha_vprop(n) > alpha_stall || M_tip_vprop(n) >= 1 || I_Bat_vprop(n) <= 0 || eta_ges_vprop(n) > 1
                C_Rest_V_vprop(n) = NaN;
                Omega_vprop(n) = NaN;
                U_mot_vprop(n) = NaN;
                I_mot_vprop(n) = NaN;
                I_Bat_vprop(n) = NaN;
                PWM_vprop(n) = NaN;
                eta_ges_vprop(n) = NaN;
                V_Kg_vprop(n) = NaN;
                
            end
            
            if x > 1     % Entferne alle Ergebnisse, wobei die Restladung im Vergleich zum vorheigen Schritt steigt
                if C_Rest_V_vprop(n) > C_Rest_V_vpp(x-1)
                    C_Rest_V_vprop(n) = NaN;
                    Omega_vprop(n) = NaN;
                    U_mot_vprop(n) = NaN;
                    I_mot_vprop(n) = NaN;
                    I_Bat_vprop(n) = NaN;
                    PWM_vprop(n) = NaN;
                    eta_ges_vprop(n) = NaN;
                    V_Kg_vprop(n) = NaN;
                end
            end
            
            
            Bestimmung_vprop(n) = Delta_C_Bat_vprop(n) * U_Bat_vprop(n);        % Berechnung der aufgebrachten Energiemenge
            
        end               % Ende der Pitchschleife
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if isnan(nanmean(Omega_vprop)) ~= 1 && isnan(nanmean(I_mot_vprop)) ~= 1 &&  isnan(nanmean(U_mot_vprop)) ~= 1 && ...     % Wenn alle der Vektoren nicht nur NaN
                isnan(nanmean(PWM_vprop)) ~= 1 && isnan(nanmean(C_Rest_V_vprop)) ~= 1 && isnan(nanmean(I_Bat_vprop)) ~= 1       % enthalten (unfliegbarer Zustand)
            
            % Kriterium für optimalen Steigwinkel
            y = 1;
            while y <= length(Bestimmung_vprop)                             % Suche und finde optimalen Steiggeschwindigkeit
                A = Bestimmung_vprop;
                A = sort(A);
                ind_opt = find(Bestimmung_vprop == A(y));
                
                if length(ind_opt) > 1
                    ind_opt = ind_opt(1);
                end

                if ~isnan(Thrust_vprop(ind_opt)) && ~isnan(Omega_vprop(ind_opt)) && ~isnan(I_mot_vprop(ind_opt)) && ~isnan(U_mot_vprop(ind_opt)) && ~isnan(PWM_vprop(ind_opt)) ...
                        && ~isnan(I_Bat_vprop(ind_opt)) && ~isnan(U_Bat_vprop(ind_opt)) && ~isnan(C_Rate_vprop(ind_opt)) && ~isnan(C_Rest_V_vprop(ind_opt))
                    
                    break;                                                 % Unterbreche Schleife, falls alle Leistungsparameter physikalisch realistische Werte besitzen
                end
                
                y = y+1;
            end
            
            % Übergabe der Leistungswerte für die optimale Geschwindigkeit und
            % Festlegen für entsprechenden Höhenschritt
            
            % Multicopter
            Theta_inter(z) = Theta_vprop(ind_opt);
            alpha_inter(z) = alpha_vprop(ind_opt);
            V_Kg_inter(z) = V_Kg_vprop(ind_opt);	    % ???????????????   Bestimmung opt. Stgeschw. und speichern in Vektor für jeden Höhenabschnitt
            % Propeller
            Thrust_inter(z) = Thrust_vprop(ind_opt);
            Omega_inter(z) = Omega_vprop(ind_opt);
            tau_inter(z) = tau_vprop(ind_opt);
            M_tip_inter(z) = M_tip_vprop(ind_opt);
            pitch_inter(z) = pitch_vprop(ind_opt);
            % Motor
            I_mot_inter(z) = I_mot_vprop(ind_opt);
            U_mot_inter(z) = U_mot_vprop(ind_opt);
            % ESC
            PWM_inter(z) = PWM_vprop(ind_opt);
            % Batterie
            I_Bat_inter(z) = I_Bat_vprop(ind_opt);
            U_Bat_inter(z) = U_Bat_vprop(ind_opt);   
            C_Rate_inter(z) = C_Rate_vprop(ind_opt);
            C_Rest_V_inter(z) = C_Rest_V_vprop(ind_opt);
            Delta_C_Bat_inter(z) = Delta_C_Bat_vprop(ind_opt);  
            i_int_inter(z) = i_int_vprop(ind_opt);  
            % Gesamtsystem
            eta_ges_inter(z) = eta_ges_vprop(ind_opt);
            
        else
            
            % Ansonsten verwerfe alle Werte
            % Multicopter
            Theta_inter(z) = NaN;
            alpha_inter(z) = NaN;
            V_Kg_inter(z) = NaN;
            % Propeller
            Thrust_inter(z) = NaN;
            Omega_inter(z) = NaN;
            tau_inter(z) = NaN;
            M_tip_inter(z) = NaN;
            pitch_inter(z) = NaN;
            % Motor
            I_mot_inter(z) = NaN;
            U_mot_inter(z) = NaN;
            % ESC
            PWM_inter(z) = NaN;
            % Batterie
            I_Bat_inter(z) = NaN;
            U_Bat_inter(z) = NaN;          
            C_Rate_inter(z) = NaN;
            C_Rest_V_inter(z) = NaN;
            % Gesamtsystem
            eta_ges_inter(z) = NaN;
            
        end
      
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Bestimmung_V_Kg(z) = Delta_C_Bat_inter(z) * U_Bat_inter(z);        % Berechnung der aufgebrachten Energiemenge
        
        z = z+1;			% Erhöhung der Zählervariablen für die gamma-Schleife
        
    end
    
    
    
    if isnan(nanmean(Omega_inter)) ~= 1 && isnan(nanmean(I_mot_inter)) ~= 1 &&  isnan(nanmean(U_mot_inter)) ~= 1 && ...     % Wenn alle der Vektoren nicht nur NaN
            isnan(nanmean(PWM_inter)) ~= 1 && isnan(nanmean(C_Rest_V_inter)) ~= 1 && isnan(nanmean(I_Bat_inter)) ~= 1       % enthalten (unfliegbarer Zustand)
        
        % Kriterium für optimalen Steigwinkel
        y = 1;
        while y < length(Bestimmung_V_Kg)+1                             % Suche und finde optimalen Steiggeschwindigkeit
            A = Bestimmung_V_Kg;
            A = sort(A);
            ind_opt2 = find(Bestimmung_V_Kg == A(y));
            
            if length(ind_opt2) > 1
                ind_opt2 = ind_opt2(1);
            end
            
            if ~isnan(Thrust_inter(ind_opt2)) && ~isnan(Omega_inter(ind_opt2)) && ~isnan(I_mot_inter(ind_opt2)) && ~isnan(U_mot_inter(ind_opt2)) && ~isnan(PWM_inter(ind_opt2)) ...
                    && ~isnan(I_Bat_inter(ind_opt2)) && ~isnan(U_Bat_inter(ind_opt2)) && ~isnan(C_Rate_inter(ind_opt2)) && ~isnan(C_Rest_V_inter(ind_opt2))
                
                break;                                                 % Unterbreche Schleife, falls alle Leistungsparameter physikalisch realistische Werte besitzen
            end
            
            y = y+1;
        end
        
        % Übergabe der Leistungswerte für die optimale Geschwindigkeit und
        % Festlegen für entsprechenden Höhenschritt
        
        % Multicopter
        Theta(x) = Theta_inter(ind_opt2);
        alpha(x) = alpha_inter(ind_opt2);
        V_Kg_vpp(x) = V_Kg_inter(ind_opt2);	    % Bestimmung opt. Stgeschw. und speichern in Vektor für jeden Höhenabschnitt
        t_Flug_vpp(x+1) = t_Flug_vpp(x) + t_Flug_inter(ind_opt2);
        % Propeller
        Thrust(x) = Thrust_inter(ind_opt2);
        Omega_vpp(x) = Omega_inter(ind_opt2);
        tau(x) = tau_inter(ind_opt2);
        M_tip(x) = M_tip_inter(ind_opt2);
        pitch_vpp(x) = pitch_inter(ind_opt2);
        % Motor
        I_mot_vpp(x) = I_mot_inter(ind_opt2);
        U_mot_vpp(x) = U_mot_inter(ind_opt2);
        % ESC
        PWM_vpp(x) = PWM_inter(ind_opt2);
        % Batterie
        I_Bat_vpp(x) = I_Bat_inter(ind_opt2);
        U_Bat_vpp(x+1) = U_Bat_inter(ind_opt2);
        C_Rate(x) = C_Rate_inter(ind_opt2);
        C_Rest_V_vpp(x) = C_Rest_V_inter(ind_opt2);
        Delta_C_Bat(x+1) = Delta_C_Bat_inter(ind_opt2);
        i_int(x+1) = i_int_inter(ind_opt2);
        % Gesamtsystem
        eta_ges_vpp(x) = eta_ges_inter(ind_opt2);
        
    else
        
        % Ansonsten verwerfe alle Werte
        % Multicopter
        Theta(x) = NaN;
        alpha(x) = NaN;
        V_Kg_vpp(x) = NaN;
        t_Flug_vpp(x+1) = NaN;
        % Propeller
        Thrust(x) = NaN;
        Omega_vpp(x) = NaN;
        tau(x) = NaN;
        M_tip(x) = NaN;
        pitch_vpp(x) = NaN;
        % Motor
        I_mot_vpp(x) = NaN;
        U_mot_vpp(x) = NaN;
        % ESC
        PWM_vpp(x) = NaN;
        % Batterie
        I_Bat_vpp(x) = NaN;
        U_Bat_vpp(x+1) = NaN;
        C_Rate(x) = NaN;
        C_Rest_V_vpp(x) = NaN;
        % Gesamtsystem
        eta_ges_vpp(x) = NaN;
        
    end
    
      
    
    H(x) = H_oben;			% Speichern der Höhe im Vektor
    x = x+1;				% Erhöhung der Zählervariablen für die Höhen-Schleife
    
    %% Spielereien
    disp([num2str((x-1)*10000/H_max) '%']);
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot

% Darstellung der Ergenisse in Diagrammen
figure(figure_ges)

subplot(621), plot(H,C_Rest_V_vpp*100,'LineWidth',2), grid, title('Restladung'), xlabel('Höhe [m]'),ylabel('C_{Bat,Rest} [%]')
subplot(622), plot(H,Omega_vpp/(2*pi)*60,'LineWidth',2), grid, title('Drehzahl'), xlabel('Höhe [m]'),ylabel('Drehzahl [RPM]')
subplot(623), plot(H,I_mot_vpp,'LineWidth',2), grid, title('Motorstrom'), xlabel('Höhe [m]'),ylabel('I_{Mot} [A]')
subplot(624), plot(H,U_mot_vpp,'LineWidth',2), grid,  title('Motorspannung'), xlabel('Höhe [m]'),ylabel('U_{mot} [V]')
subplot(625), plot(H,I_Bat_vpp,'LineWidth',2), grid, title('Batteriestrom'), xlabel('Höhe [m]'),ylabel('I_{Bat} [A]')
H2 = [0;H];
subplot(626), plot(H2,U_Bat_vpp,'LineWidth',2), grid, title('Batteriespannung'), xlabel('Höhe [m]'),ylabel('U_{Bat} [V]')
subplot(627), plot(H,PWM_vpp*100,'LineWidth',2), grid, title('Pulsweitenmodulation'), xlabel('Höhe [m]'),ylabel('PWM [%]')
subplot(628), plot(H,eta_ges_vpp*100,'LineWidth',2), grid, title('Gesamtwirkungsgrad'), xlabel('Höhe [m]'),ylabel('eta_{ges} [%]')
subplot(629), plot(H,V_Kg_vpp,'LineWidth',2), title('Bahngeschwindigkeit'), grid, xlabel('Höhe [m]'),ylabel('V_{Kg} [m/s]')
subplot(6,2,10), plot(H2,t_Flug_vpp,'LineWidth',2), title('Flugzeit'), grid, xlabel('Höhe [m]'),ylabel('t_{Flug} [s]')
subplot(6,2,11), plot(H,pitch_vpp,'LineWidth',2), title('Pitch'), grid, xlabel('Höhe [m]'),ylabel('Pitch [in]')
% Anpassung und Abspeichern der Diagramme
PaperSizeX = 21;
PaperSizeY = 29.7;

fig = gcf;
set(gcf,'PaperUnits','centimeters', 'PaperPosition', [-1.75 -2.1 24.65 34.45]);%[-1.75 -2.6 24.65 34.45]);
set(gcf,'Units','centimeters', 'PaperSize', [PaperSizeX PaperSizeY]);
saveas(gcf,Dateiname, 'pdf');


