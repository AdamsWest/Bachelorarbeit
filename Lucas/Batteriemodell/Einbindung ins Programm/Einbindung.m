% Ablauf in Leistungsberechnung

%% Intialisierungen
  
Q = 1.0618;                   % Q
Qnom = 0.9102;                % Qnom
Qexp = 0.2083;                % Qexp
Vfull = 4.0399;               % Vfull
Vexp = 3.7783;                % Vexp
Vnom = 3.5181;                % Vnom
i = 29.7754;                  % i
R_bat = 0.0142;               % R_bat
Batterie_data = [Q Qnom Qexp Vfull Vexp Vnom i R_bat 0 0 0];        % Zwischenbelegung
[Eo,A,K] = Batterie_parameter(Batterie_data);
Batterie_data = [Q Qnom Qexp Vfull Vexp Vnom i R_bat Eo A K];       % vollständiger Vektor


% Initialisierungen
i_int = zeros(lengthi,1);

%% In der for-Schleife, Funktionsaufruf

[I_bat,U_bat,C_rate,Delta_C_bat,C_Rest_V,i_int(i)] = Batterie(PWM,eta_PWM,I_mot,n_Prop,C_bat,P_bat,Delta_C_bat,t_Flug,i_int(i-1),N_Bat_cell,Batterie_data);


%% Funktion

% function [I_bat,U_bat,C_rate,Delta_C_bat,C_Rest_V,i_int] = Batterie(PWM,eta_PWM,I_mot,n_Prop,C_bat,P_bat,Delta_C_bat,t_Flug,i_int_ante,N_Bat_cell,Batterie_data,)

% B = 3/Batterie_data(3);
% Q = Batterie_data(1);
% Qexp = Batterie_data(3);
% R = Batterie_data(8);
% Eo = Batterie_data(9);
% A = Batterie_data(10);
% K = Batterie_data(11);
% i_star = 0;

% I_bat = PWM * I_mot / eta_PWM * n_Prop;                     % Batteriestrom
% i_int = I_bat*t_Flug + i_int_ante;                               % Integral des Stroms
% U_bat = Eo - R*I_bat - K * Q / (Q - i_int) * (i_int_post + I_bat*i_star) ...
%    + A * exp(-B*i_int);                                     % Batteriespannung pro Zelle
% U_bat = N_Bat_cell * U_bat;                                 % gesamt Batteriespannung

% C_rate = I_bat / (C_bat/3600);                              % C-Rate bezogen auf eine nominale Entladezeit von 1 Stunde
% C_bat_Peukert = C_bat * (1/C_rate)^(P_bat-1);               % nutzbare Kapazitaet nach Peukert
% Delta_C_bat = I_bat * t_Flug + Delta_C_bat;                 % entnommene Ladung
% C_Rest_V = (C_bat_Peukert - Delta_C_bat) / C_bat_Peukert;   % Ladezustand als Verhältnis
%%



