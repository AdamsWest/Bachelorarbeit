% Leistungsberechnung für Flugsysteme

%% Initialisierungen

% Batterie

U_Bat_nom = N_Bat_cell * U_Bat_cell;        % nominale Batteriespannung
U_Bat_min = N_Bat_cell * U_Bat_cell_min;    % minimale Batteriespannung
C_Bat = E_Dichte * m_Bat / U_Bat_nom;       % Kapazitaet der Batterie in As
Delta_C_Bat = 0;                            % Initialisierung Batteriekapazität, die nach jedem delta_h gebraucht wird



% Propeller

R = D * 0.0254 / 2;                         % Propellerradius in Meter
F = pi * R^2;                               % Fläche eines Propellers in Quadratmeter
Theta_75 = atan( 4*P_75 / (3*pi * D) );     % geometrischer Anstellwinkel des Propellers bei 75% des Radius
[RPM_map, V_map, T_map, P_map, TAU_map] = Propeller_map(DATA_APC,prop_name);    % Aufbau des Kennfeldes



% Copter

m = m_copter; %+ m_Bat;                     % Gesamtmasse des Quadrocopters



% Umgebungsparameter

T_11 = T_0 - 0.0065 *(11000-H_0);                   % T in 11000m Höhe                     
rho_11 = rho_0 * (1 - 0.0065*(11000/T_0))^4.256;    % Dichte in 11000m Höhe
u_Wg = V_Wg(1);                             % Windgeschwindigkeiten definieren
v_Wg = V_Wg(2);
w_Wg = V_Wg(3);



% Intialisierung der Matrizen für jeden Höhenabschnitt

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
%v = zeros(lengthi,1);

t_Flug = Delta_H / V_Kg;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Programmanfang

i = 1;                                      % Zähler intialisiern
for h = H_0:Delta_H:H_max
    
    % Berechnung der Flughöhe für iterative Schritte und der mittleren
    % Dichte zwischen den Diskretisierungspunkten
    
    H_unten = h;
    H_oben = h + Delta_H;
    
    
    % Berechnung der Dichte an den Intevallgrenzen
    
    if H_oben <= 11000
        rho_unten = rho_0 *(1-0.0065*(H_unten/T_0))^4.256;
        rho_oben = rho_0 *(1-0.0065*(H_oben/T_0))^4.256;
    elseif H_unten <= 11000 && H_oben >11000
        rho_unten = rho_0 *(1 - 0.0065*(H_unten/T_0))^4.256;
        rho_oben = rho_11 * exp(-g/(287*T_11)*(H_oben-11000));
    else
        rho_unten = rho_11 * exp(-g/(287*T_11)*(H_unten-11000));
        rho_oben = rho_11 * exp(-g/(287*T_11)*(H_oben-11000));
    end
    
    
    rho(i) = rho_unten + (rho_oben - rho_unten)/2;                  % Berechnung der mittleren Dichte im Intervall
    
    
    if i == 1
        rho_1 = rho_0;
        rho_2 = rho(i);
    else
        rho_1 = rho(i-1);
        rho_2 = rho(i);
    end
    
    
    T_map = T_map * rho(i)/rho_1;                                   % Anpassung des Schubkennfeldes an die sich ändernde Dichte
    P_map = P_map * rho(i)/rho_1;                                   % Anpassung des Leistungskennfeldes an die sich ändernde Dichte
    TAU_map = TAU_map * rho(i)/rho_1;
    
    
    %switch Vorgabe
        %case 'V'
            
            if Abfrage_Flugsystem == 1                                      % handelt es sich um Multicopter oder Flächenflugzeug
                
                % MULTICOPTER
                % Aerodynamik berechnen
                
                [Thrust(i),Theta(i),V_A,alpha(i)] = MulticopterAerodynamik(u_Wg,V_Kg,gamma,c_W_copter_seitlich,c_W_copter_oben,c_A_copter_max,rho(i),A_copter,m,g);
                
            else
                
                % FLÄCHENFLUGZEUG
                
                % Aerodynamik
                
                [Thrust] = FlaechenflugzeugAerodynamik(m,g,epsilon);
                
                
            end
            
            
            Thrust(i) = Thrust(i) / n_Prop;                                 % Schub auf n Propeller verteilen
            
            
            
            if Thrust(i) > max(max(T_map))                                  % wenn, Schub zu gross (Ergebnis verwerfen)
                Omega(i) = NaN;
                I_mot(i) = NaN;
                C_Rate(i) = NaN;
                C_Rest_V(i) = NaN;
            else
                
                
                % Drehzahl und Drehmoment bestimmen
                
                [Omega(i),tau(i)] = Propeller(V_A, alpha(i), Thrust(i), RPM_map, V_map, T_map, TAU_map);
                
                % Motorzustand berechnen
                
                [U_mot(i),I_mot(i)] = Motor(tau(i),K_V,I_0,R_i,Omega(i));
                
                
                % Zustand der Motorregler berechnen
                
                [PWM,eta_PWM] = ESC(U_mot(i),U_Bat_nom);
                
                
                % Batteriezustand berechnen
                
                [I_Bat,C_Rate(i),Delta_C_Bat,C_Rest_V(i)] = Batterie(PWM,eta_PWM,I_mot(i),n_Prop,C_Bat,P_Bat,Delta_C_Bat,t_Flug);
            end
        %case 'PWM'
            
            % hier Flugprofil einfügen und Funktionen umschreiben
            %if H_unten <= 1300
            %    PWM = PWM_Profil(1);
            %elseif H_unten > 1300 && H_unten <= 6700
            %    PWM = PWM_Profil(2);
            %elseif H_unten > 6700 && H_unten <= 8800
            %    PWM = PWM_Profil(3);
            %elseif H_unten > 8800 && H_unten <= 10000
            %    PWM = PWM_Profil(4);
            %else
            %    PWM = PWM_Profil(5);
            %end
            
            
            % Motorstrom berechnen
            
            [U_mot(i),eta_PWM] = ESC_umgk(PWM,U_Bat_nom);
            
            
            [I_mot(i),Omega(i),tau(i)] = Motor_umgk(PWM,U_mot(i),P_map,I_0,K_V,R_i);
            
            
            [v(i)] = Propeller_umgk(Omega(i),PWM,V_map,RPM_map,T_map);
            
            
            [I_Bat,C_Rate(i),Delta_C_Bat,C_Rest_V(i)] = Batterie(PWM,eta_PWM,I_mot(i),n_Prop,C_Bat,P_Bat,Delta_C_Bat,t_Flug);
            
    end
    
    
    
    % Werden Grenzen ueberschritten?
    
    
    % Wenn Grenzen ueberschritten werden, Resultate entfernen
    
    if C_Rest_V(i) < 0.0 || U_mot(i) > U_Bat_nom || U_mot(i) <= 0 || C_Rate(i) > C_Rate_max || I_mot(i) > I_max || alpha(i) > alpha_stall
        C_Rest_V(i) = NaN;
        Omega(i) = NaN;
        U_mot(i) = NaN;
        I_mot(i) = NaN;
    end
    
    
    H(i) = H_oben;
    I_bat_ges(i) = I_Bat;
    i = i+1;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot

% Restladung über der Höhe
figure(figure_C_Rest_V)
plot(H,C_Rest_V*100,'LineWidth',2);
grid on
hold on
xlabel('Höhe [m]')
ylabel('Restladung [%]')


% Drehzahl über der Höhe
figure(figure_omega)
plot(H,Omega/(2*pi)*60,'LineWidth',2)
grid on
hold on
xlabel('Höhe [m]')
ylabel('RPM')

% Motorstrom über der Höhe
figure(figure_I_mot)
plot(H,I_mot,'LineWidth',2)
grid on
hold on
xlabel('Höhe [m]')
ylabel('I_{mot} [A]')

% Motorspannung über der Höhe
figure(figure_U_mot)
plot(H,U_mot,'LineWidth',2)
grid on
hold on
xlabel('Höhe [m]')
ylabel('U_{mot} [V]')

% Batteriestrom über der Höhe
figure(figure_I_Bat)
plot(H,I_bat_ges,'LineWidth',2)
grid on
hold on
xlabel('Höhe [m]')
ylabel('I_{Bat} [V]')
