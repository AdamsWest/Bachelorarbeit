function [I_bat,U_bat,C_rate,Delta_C_bat,C_Rest_V,i_int] = Batterie_neu(Batterie_data,Cnom,PWM,eta_PWM,n_Prop,i_int,U_bat,C_bat,Delta_C_bat,I_mot,N_bat_cell,P_bat,t_Flug)


% BATTERIE   Berechnet die Restladung der Batterie in Abhängigkeit des Motorstroms
%
%   Batterie berchnet anhand der PWM Stellung, des Motorstroms, der
%   Flugzeit und anderer Parameter die restliche Ladung, den Batteriestrom
%   und die Entladerate
%
% Syntax:  [I_bat,C_rate,Delta_C_bat,C_Rest_V] = Batterie(PWM,eta_PWM,I_mot,n_Prop,C_bat,P_bat,Delta_C_bat,t_Flug)
%
% Inputs:
%   PWM     Pulsweitenmodulation (Verhältnis Motorspannung zur maximal
%   verfügbaren Spannung
%   eta_PWM Wirkungsgrad des Motorreglers
%   I_mot   Motorstrom
%   n_Prop  Propelleranzahl
%   C_bat   Batteriekapazität
%   P_bat   Peukert-Konstante
%   Delta_C_bat     entnommene Ladung
%   t_Flug  Flugzeit
%
% Outputs:
%   I_bat   Batteriestrom
%   C_rate  C-Rate
%   Delta_C_bat     neu berechenete entnommene Ladung
%   C_Rest_V    verbleibende Kapazität
%
% Example:
%    Line 1 of example
%    Line 2 of example
%    [ values, derivatives ] = Untitled( x, myStruct.y, 2 )
%    [ values, ~ ] = Untitled( [ 1:.1:100 ] , myStruct.y, 1 )
%
% See also: ESC

%   Copyright 2018 TU-Braunschweig
% ******************************************************************************

% Übergeben wird: Batterie_data, Cnom, i_int, PWM, eta_PWM, I_mot, n_Prop,
% N_bat_cell, 
% Was wird ausgegeben: I_bat, U_bat, C_rate, Delta_C_bat, C_Rest_V, i_int



step_dt = 1;                                                % Zeitschritt Delta_t in Sekunden festlegen 

% hier wird Vektor mit Daten übergeben
% Batterie_data = [Q Qnom Qexp Vfull Vexp Vnom i 0 0 0];        % Zwischenbelegung
Q = Batterie_data(1);
Qnom = Batterie_data(2);
Qexp = Batterie_data(3);
Vfull = Batterie_data(4);
Vexp = Batterie_data(5);
Vnom = Batterie_data(6);
i = Batterie_data(7);
B = 3/Qexp;


% Berechnung weiterer Parameter für die Normzelle

I_bat = PWM * I_mot / eta_PWM * n_Prop;                     % Batteriestrom
I_bat = I_bat/Cnom;                                         % Normierung
C_rate = I_bat / (C_bat/3600);                              % C-Rate bezogen auf eine nominale Entladezeit von 1 Stunde
R = 0.1077 / (Q*(0.1555*Q+0.9825*C_rate)^0.5485);           % Modell für den Innenwiderstand, Parameter nach geringster Fläche ausgewählt

Batterie_data = [Q Qnom Qexp Vfull Vexp Vnom i R];        % Zwischenbelegung
[Eo,A,K] = Batterie_parameter(Batterie_data);

% i_int = 0;
% U_bat = 0;


for n = 1:floor(t_Flug)                                            % Iteration für über der Zeit, mit der mit der konstanten C_Rate geflogen wird
    
    
    i_int = I_bat*step_dt/3600 + i_int;   % integral of the current
    
    
    % calculating the battery voltage (of one cell)
    U_bat = Eo - R*I_bat - K * Q / (Q - i_int) * (i_int + I_bat*0) + A * exp(-B*i_int);
    % the battery voltage of all cells
    U_bat = N_bat_cell * U_bat;
    
    
end

% Berechnung der Restkapazität

I_bat = PWM * I_mot / eta_PWM * n_Prop;                     % Batteriestrom (Neu- bzw. Wiederbelegung)
C_rate = I_bat / (C_bat/3600);                              % C-Rate bezogen auf eine nominale Entladezeit von 1 Stunde (Neu- bzw. Wiederbelegung)
C_bat_Peukert = C_bat * (1/C_rate)^(P_bat-1);               % nutzbare Kapazitaet nach Peukert
Delta_C_bat = I_bat * t_Flug + Delta_C_bat;                 % entnommene Ladung
C_Rest_V = (C_bat_Peukert - Delta_C_bat) / C_bat_Peukert;   % Ladezustand als Verhältnis


end

