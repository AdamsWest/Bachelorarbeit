%% Leistungsberechnung f�r Flugsysteme

%% Initialisierungen

% Batterie

U_Bat_nom = N_Bat_cell * U_Bat_cell;        % nominale Batteriespannung
U_Bat_min = N_Bat_cell * U_Bat_cell_min;    % minimale Batteriespannung
C_Bat = E_Dichte * m_Bat / U_Bat_nom;       % Kapazitaet der Batterie in As
Delta_C_Bat = 0;                            % Initialisierung Batteriekapazit�t, die nach jedem delta_h gebraucht wird



% Propeller

  
R = D * 0.0254 / 2;                         % Propellerradius in Meter
F = pi * R^2;                               % Fl�che eines Propellers in Quadratmeter
%Theta_75 = atan( 4*P_75 / (3*pi * D) );     % geometrischer Anstellwinkel des Propellers bei 75% des Radius
%[RPM_map, V_map, T_map, P_map, TAU_map] = Propeller_map(DATA_APC,prop_name);    % Aufbau des Kennfeldes

DATA = evalin('base','DATA_APC');

[l,w] = size(DATA_APC);
%% 
% Aufsuchen und finden aller Propellernamen, die den gleichen Durchmesser
% wie der gesuchte haben. Dabei werden alle anderen Propeller entfernt

for m = l:-1:1
    a = DATA{m,1};                                  % Speichern des Propellernames unter der Variablen a
    diameter = str2double(a(1:strfind(a,'x')-1));   % Extrahieren des Durchmessers aus dem Propellernamen
    if diameter ~= D                                % Entfernen der Zeile mit einem anderen Propellerdurchmesser
        DATA(m,:) = [];
    end
end

%% 
% Entnahme jedes einzelnen Propellers
counter = 0;
ind = 1;
len = length(DATA);
d = 1;
while ind < len
    
    prop_name = DATA{ind,1};
    pitch_all = prop_name(strfind(prop_name,'x')+1:end);
    ind = find(strcmp(DATA(:,1),prop_name));
    
    
    ind_unten = min(ind);
    ind_oben = max(ind);
    ind = max(ind) + 1;
    
    
    % nicht passende Propeller werden entfernt (-E, 3- and 4-blade, -SF, etc.)
    if isnan(str2double(pitch_all)) == 1
        
        for d = ind_oben:-1:ind_unten
            DATA(d,:) = [];
            len = len - 1;
        end
        
    else
        
        [RPM, V, T, P, Tau] = Propeller_map(DATA,prop_name);
        
        
        counter = counter + 1;
        
        % Die Daten f�r RPM, V, T, usw. werden fortlaufenden Vektoren
        % zugewiesen
        assignin ('base',['prop_name' num2str(counter)], prop_name);
        assignin ('base',['pitch_' num2str(counter)], pitch_all);
        assignin ('base',['RPM_' num2str(counter)], RPM);
        assignin ('base',['T_' num2str(counter)], T);
        assignin ('base',['P_' num2str(counter)], P);
        assignin ('base',['Tau_' num2str(counter)], Tau);
    end
end



% Umgebungsparameter

Temp_11 = Temp_0 - 0.0065 *(11000-H_0);                   % T in 11000m H�he                     
rho_11 = rho_0 * (1 - 0.0065*(11000/Temp_0))^4.256;    % Dichte in 11000m H�he




% Intialisierung der Matrizen f�r jeden H�henabschnitt

lengthi = floor(abs(H_max - H_0) / Delta_H + 1); 
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
I_bat_ges = zeros(lengthi,1);
PWM = zeros(lengthi,1);
Pitch = zeros(lengthi,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Programmanfang

i = 1;                                      % Z�hler intialisiern
for h = H_0:Delta_H:H_max
    
    % Berechnung der Flugh�he f�r iterative Schritte und der mittleren
    % Dichte zwischen den Diskretisierungspunkten
    
    H_unten = h;
    H_oben = h + Delta_H;
    
    
    % Berechnung der Dichte an den Intevallgrenzen
    
    if H_oben <= 11000
        
        rho_unten = rho_0 *(1-0.0065*(H_unten/Temp_0))^4.256;
        rho_oben = rho_0 *(1-0.0065*(H_oben/Temp_0))^4.256;
        
    elseif H_unten <= 11000 && H_oben >11000
        
        rho_unten = rho_0 *(1 - 0.0065*(H_unten/Temp_0))^4.256;
        rho_oben = rho_11 * exp(-g/(287*Temp_11)*(H_oben-11000));
        
    else
        
        rho_unten = rho_11 * exp(-g/(287*Temp_11)*(H_unten-11000));
        rho_oben = rho_11 * exp(-g/(287*Temp_11)*(H_oben-11000));
        
    end
    
    
    rho(i) = rho_unten + (rho_oben - rho_unten)/2;                  % Berechnung der mittleren Dichte im Intervall
    
    
    if i == 1
        rho_1 = rho_0;
        rho_2 = rho(i);
    else
        rho_1 = rho(i-1);
        rho_2 = rho(i);
    end
    
    
    %T_map = T_map * rho(i)/rho_1;                                   % Anpassung des Schubkennfeldes an die sich �ndernde Dichte
    %P_map = P_map * rho(i)/rho_1;                                   % Anpassung des Leistungskennfeldes an die sich �ndernde Dichte
    %TAU_map = TAU_map * rho(i)/rho_1;                               % Anpassung des Drehmomentkennfeldes an die sich �ndernde Dichte
    
    
    % Steiggeschwindigkeitsprofil vorgeben
    %if H_unten < 300
    %    V_Kg = V_Profil(1);
    %elseif H_unten >= 300 && H_unten < 1700
    %    V_Kg = V_Profil(2);
    %elseif H_unten >= 1700 && H_unten < 3100
    %    V_Kg = V_Profil(3);
    %elseif H_unten >= 3100 && H_unten < 5000
    %    V_Kg = V_Profil(4);
    %elseif H_unten >= 5000 && H_unten < 6000
    %    V_Kg = V_Profil(5);
    %elseif H_unten >= 6000 && H_unten < 7700
    %    V_Kg = V_Profil(6);
    %elseif H_unten >= 7700 && H_unten < 9800
    %    V_Kg = V_Profil(7);
    %elseif H_unten >= 9800 && H_unten < 10300
    %    V_Kg = V_Profil(8);
    %else
    %    V_Kg = 3;
    %end
       
    t_Flug = Delta_H / V_Kg;                                        % Fluggeschwindigkeit
        
    if Abfrage_Flugsystem == 1                                      % handelt es sich um Multicopter oder Fl�chenflugzeug
        
        % MULTICOPTER
        m = m_copter + m_Bat + m_Mot * n_Prop;                     % Gesamtmasse des Quadrocopters
        
        % Aerodynamik berechnen
        
        [Thrust(i),Theta(i),V_A,alpha(i)] = MulticopterAerodynamik(u_Wg,V_Kg,gamma,c_W_copter_seitlich,c_W_copter_oben,c_A_copter_max,rho(i),A_copter,m,g);
        
    else
        
        % FL�CHENFLUGZEUG
        m = m_flugzeug + m_Bat + m_Mot * n_Prop;                  % Gesamtmasse des Fl�chenflugzeugs
        
        % Aerodynamik
        
        [Thrust(i),V_A] = FlaechenflugzeugAerodynamik(m,g,epsilon,V_Kg);
        
        
    end
    
    
    Thrust(i) = Thrust(i) / n_Prop;                                 % Schub auf n Propeller verteilen
    
    
    
    if Thrust(i) > max(max(T_1))            % !!!!!!!!!!!!! Ab�ndern                      % wenn Schub zu gross (Ergebnis verwerfen)
        Omega(i) = NaN;
        I_mot(i) = NaN;
        C_Rate(i) = NaN;
        C_Rest_V(i) = NaN;
        I_bat_ges(i) = NaN;
        PWM(i) = NaN;
        Pitch(i) = NaN;
    else
        
        
        % Drehzahl und Drehmoment bestimmen
        
        %[Omega(i),tau(i)] = Propeller(V_A, alpha(i), Thrust(i), RPM_map, V_map, T_map, TAU_map);
        
        omega = zeros(counter,1);
        TAU = zeros(counter,1);
        
        %% Bestimmung des Pitches mit der geringsten Leistung
        % Theoretisch in gleicher While-Schleife m�glich wie in
        pitch_all = zeros(counter,1);
        for n = 1:counter
            RPM = evalin('base',['RPM_' num2str(n)]);
            T = evalin('base',['T_' num2str(n)]) * rho_2/rho_1;
            Tau = evalin('base',['Tau_' num2str(n)]) * rho_2/rho_1;
            
            a = str2double(evalin('base',['pitch_' num2str(n)]));
            pitch_all(n) = a;
            
            [omega(n),TAU(n)] = Propeller(V_A, alpha(n), Thrust(n), RPM, V, T, Tau);
        end
        
        tau(i) = min(TAU);
        mintorqueind = find(TAU == tau(i));
        Omega(i) = omega(mintorqueind);
        Pitch(i) = pitch_all(mintorqueind);
        %mintorquepitch = evalin('base',['pitch_' num2str(mintorqueind)]);
        
        
        
        % Motorzustand berechnen
        
        [U_mot(i),I_mot(i)] = Motor(tau(i),K_V,I_0,R_i,Omega(i));
        
        
        % Zustand der Motorregler berechnen
        
        [PWM(i),eta_PWM] = ESC(U_mot(i),U_Bat_nom);
        
        
        % Batteriezustand berechnen
        
        [I_Bat,C_Rate(i),Delta_C_Bat,C_Rest_V(i)] = Batterie(PWM(i),eta_PWM,I_mot(i),n_Prop,C_Bat,P_Bat,Delta_C_Bat,t_Flug);
        I_bat_ges(i) = I_Bat;
        
        
    end
    
    
    
    
    % Werden Grenzen ueberschritten?
    
    
    % Wenn Grenzen ueberschritten werden, Resultate entfernen
    
    if C_Rest_V(i) < 0.0 || U_mot(i) > U_Bat_nom || U_mot(i) <= 0 || C_Rate(i) > C_Rate_max || I_mot(i) > I_max || alpha(i) > alpha_stall
        C_Rest_V(i) = NaN;
        Omega(i) = NaN;
        U_mot(i) = NaN;
        I_mot(i) = NaN;
        I_bat_ges(i) = NaN;
        PWM(i) = NaN;
        Pitch(i) = NaN;
    end
    
    
    H(i) = H_oben;
    i = i+1;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot

% Restladung �ber der H�he
figure(figure_C_Rest_V)
plot(H,C_Rest_V*100,'LineWidth',2);
grid on
hold on
xlabel('H�he [m]')
ylabel('Restladung [%]')


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
plot(H,I_bat_ges,'LineWidth',2)
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

% Verwendeter Pitch
figure(figure_Pitch)
plot(H,Pitch,'rx')
grid on
hold on
xlabel('H�he [m]')
ylabel('Pitch [in]')

%% Datei abspeichern
%ImageSizeX = 14;
%ImageSizeY = 24;
%figure(figure_C_Rest_V)
%set(gcf,'PaperUnits','centimeters', 'PaperPosition', [0 0 ImageSizeX ImageSizeY]); 
%set(gcf,'Units','centimeters', 'PaperSize', [ImageSizeX ImageSizeY]); 
%saveas(gcf,Dateiname, 'pdf'); 