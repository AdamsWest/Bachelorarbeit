%% Leistungsberechnung f�r Flugsysteme %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% *************************************************************************


%% Initialisierungen %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Batterie

U_Bat_nom = N_Bat_cell * U_Bat_cell;        % nominale Batteriespannung
U_Bat_min = N_Bat_cell * U_Bat_cell_min;    % minimale Batteriespannung
% C_Bat = N_Bat_cell_p*C_Bat_cell*3600;       % Kapazit�t in As
% Delta_C_Bat = 0;                            % Initialisierung Batteriekapazit�t, die nach jedem delta_h gebraucht wird
C_Bat = E_Dichte * m_Bat / U_Bat_nom;                               % Kapazitaet der Batterie in As

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

% Matrixl�ngen
lengthi = floor(abs(H_max - H_0) / Delta_H + 1);
lengthj = length(Abfrage_getriebe);
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


for j = 1:lengthj
    
    m_getriebe = Abfrage_getriebe(j);
    
    m = m_copter + m_Bat + (m_Mot + m_getriebe) * n_Prop + m_nutz;                     % Gesamtmasse des Quadrocopters
    


    
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
        
        % Gesamtl�nge der Vektoren: lengthall = lengthgamma + lengthb_vert;
        
        
        
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
        for V_Kg_variabel = V_Kg_min:V_Kg_Delta:V_Kg_max	    % Variation des Bahnneigungswinkel f�r das Fl�chenflugzeug
            
            V_Kg_inter(z) = V_Kg_variabel;
            
            % MULTICOPTER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                       
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
                
                
                % Wie gro� ist die Blattspitzengeschwindigkeit?
                M_tip_inter(z) = (Omega_inter(z) * R)/a;                       % Blattspitzengeschwindigkeit in Ma
                
                
                % Motorzustand berechnen
                [U_mot_inter(z),I_mot_inter(z)] = Motor(tau_inter(z),K_V,I_0,R_i,Omega_inter(z));
                
                
                % Zustand der Motorregler berechnen
                [PWM_inter(z),eta_PWM] = ESC(U_mot_inter(z),U_Bat(x));         % <-- hier U_bat_inter
                
                
                % Batteriezustand berechnen
                Delta_C_Bat_inter(z) = Delta_C_Bat(x);                         % Anpassung der Batteriekapazit�t
                i_int_inter(z) = i_int(x);                                     % �bergabe des Integrals der Spannung vom letzten Schritt
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
                
                eta_ges_inter(z) = (n_Prop * Thrust_inter(z) * (mu_z + vi))/(I_Bat_inter(z) * U_Bat_nom);         % Leistung, die in Schub umgesetzt wird im Verh�ltnis zur aufgebrachten Leistung
                
                
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
                V_Kg_inter(z) = NaN;
                
            end
            
            
            Bestimmung_V_Kg(z) = Delta_C_Bat_inter(z) * U_Bat_inter(z);        % Berechnung der aufgebrachten Energiemenge
            
            z = z+1;			% Erh�hung der Z�hlervariablen f�r die gamma-Schleife
            
        end
        
        
        if isnan(nanmean(Omega_inter)) ~= 1 && isnan(nanmean(I_mot_inter)) ~= 1 &&  isnan(nanmean(U_mot_inter)) ~= 1 && ...     % Wenn alle der Vektoren nicht nur NaN
                isnan(nanmean(PWM_inter)) ~= 1 && isnan(nanmean(C_Rest_V_inter)) ~= 1 && isnan(nanmean(I_Bat_inter)) ~= 1       % enthalten (unfliegbarer Zustand)
            
            % Kriterium f�r optimalen Steigwinkel
            y = 1;
            while y < length(Bestimmung_V_Kg)                           % Suche und finde optimalen Steiggeschwindigkeit
                A = Bestimmung_V_Kg;                                    
                A = sort(A);                                            % Belegung A mit dem sortierten Werten
                ind_opt = find(Bestimmung_V_Kg == A(y));                % Index mit optimalen Flugzustand
                
                if length(ind_opt) > 1                                  % Wenn die L�nge des Index gr��er als 1 ist ...
                    ind_opt = ind_opt(1);                               % nimm den ersten Eintrag
                end
                
                if ~isnan(Thrust_inter(ind_opt)) && ~isnan(Omega_inter(ind_opt)) && ~isnan(I_mot_inter(ind_opt)) && ~isnan(U_mot_inter(ind_opt)) && ~isnan(PWM_inter(ind_opt)) ...
                        && ~isnan(I_Bat_inter(ind_opt)) && ~isnan(U_Bat_inter(ind_opt)) && ~isnan(C_Rate_inter(ind_opt)) && ~isnan(C_Rest_V_inter(ind_opt))
                    % Falls alle Eintr�ge der Vektoren mit ind_opt real
                    % sind, dann ...
                    break;                                                 % Unterbreche Schleife, falls alle Leistungsparameter physikalisch realistische Werte besitzen
                end
                
                y = y+1;
            end
            
            % �bergabe der Leistungswerte f�r die optimale Geschwindigkeit und
            % Festlegen f�r entsprechenden H�henschritt
            
            % Multicopter
            Theta(x) = Theta_inter(ind_opt);
            alpha(x) = alpha_inter(ind_opt);
            V_Kg(x) = V_Kg_inter(ind_opt);	    % Bestimmung opt. Stgeschw. und speichern in Vektor f�r jeden H�henabschnitt
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
        
        
        H(x) = H_oben;			% Speichern der H�he im Vektor
        x = x+1;				% Erh�hung der Z�hlervariablen f�r die H�hen-Schleife
        
    end
    
    % Darstellung der Ergenisse in Diagrammen
    figure(figure_ges)
    
    subplot(521), l1(j) = plot(H,C_Rest_V*100,'LineWidth',2); l1_Info{j} = ['m_{Getriebe} = ' num2str(m_getriebe) ' kg']; grid on, hold on
    subplot(522), l2(j) = plot(H,Omega/(2*pi)*60,'LineWidth',2); l2_Info{j} = ['m_{Getriebe} = ' num2str(m_getriebe) ' kg']; grid on, hold on
    subplot(523), l3(j) = plot(H,I_mot,'LineWidth',2); l3_Info{j} = ['m_{Getriebe} = ' num2str(m_getriebe) ' kg']; grid on, hold on
    subplot(524), l4(j) = plot(H,U_mot,'LineWidth',2); l4_Info{j} = ['m_{Getriebe} = ' num2str(m_getriebe) ' kg']; grid on, hold on
    subplot(525), l5(j) = plot(H,I_Bat,'LineWidth',2); l5_Info{j} = ['m_{Getriebe} = ' num2str(m_getriebe) ' kg']; grid on, hold on
    H2 = [0;H];
    subplot(526), l6(j) = plot(H2,U_Bat,'LineWidth',2); l6_Info{j} = ['m_{Getriebe} = ' num2str(m_getriebe) ' kg']; grid on, hold on
    subplot(527), l7(j) = plot(H,PWM*100,'LineWidth',2); l7_Info{j} = ['m_{Getriebe} = ' num2str(m_getriebe) ' kg']; grid on, hold on
    subplot(528), l8(j) = plot(H,eta_ges*100,'LineWidth',2); l8_Info{j} = ['m_{Getriebe} = ' num2str(m_getriebe) ' kg']; grid on, hold on
    subplot(529), l9(j) = plot(H,V_Kg,'LineWidth',2); l9_Info{j} = ['m_{Getriebe} = ' num2str(m_getriebe) ' kg']; grid on, hold on
    
    %% Spielereien
    disp([num2str((j)*100/lengthj) '%']);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot

% Darstellung der Ergenisse in Diagrammen
figure(figure_ges)

subplot(521), title('Restladung'), xlabel('H�he [m]'),
    ylabel('C_{Bat,Rest} [%]'), legend([l1], l1_Info)
subplot(522), title('Drehzahl'), xlabel('H�he [m]'),
    ylabel('Drehzahl [RPM]'), legend([l2], l2_Info)
subplot(523), title('Motorstrom'), xlabel('H�he [m]'),
    ylabel('I_{Mot} [A]'), legend([l3], l3_Info)
subplot(524), title('Motorspannung'), xlabel('H�he [m]'),
    ylabel('U_{mot} [V]'), legend([l4], l4_Info)
subplot(525), title('Batteriestrom'), xlabel('H�he [m]'),
    ylabel('I_{Bat} [%]'), legend([l5], l5_Info)
H2 = [0;H];
subplot(526), title('Batteriespannung'), xlabel('H�he [m]'),
    ylabel('U_{Bat} [A]'), legend([l6], l6_Info)
subplot(527), title('Pulsweitenmodulation'), xlabel('H�he [m]'),
    ylabel('PWM [%]'), legend([l7], l7_Info)
subplot(528), title('Gesamtwirkungsgrad'), xlabel('H�he [m]'),
    ylabel('eta_{ges} [%]'), legend([l8], l8_Info)
subplot(529), title('Bahngeschwindigkeit'), xlabel('H�he [m]'),
    ylabel('V_{Kg} [m/s]'), legend([l9], l9_Info)

% Anpassung und Abspeichern der Diagramme
ImageSizeX = 14;
ImageSizeY = 24;
figure(figure_ges)
set(gcf,'PaperUnits','centimeters', 'PaperPosition', [0 0 ImageSizeX ImageSizeY]);
set(gcf,'Units','centimeters', 'PaperSize', [ImageSizeX ImageSizeY]);
saveas(gcf,Dateiname, 'pdf');




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