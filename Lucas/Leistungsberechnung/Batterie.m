function [I_bat,C_rate,Delta_C_bat,C_Rest_V] = Batterie(PWM,eta_PWM,I_mot,n_Prop,C_bat,P_bat,Delta_C_bat,t_Flug)

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

I_bat = PWM * I_mot / eta_PWM * n_Prop;                     % Batteriestrom
C_rate = I_bat / (C_bat/3600);                              % C-Rate bezogen auf eine nominale Entladezeit von 1 Stunde
C_bat_Peukert = C_bat * (1/C_rate)^(P_bat-1);               % nutzbare Kapazitaet nach Peukert
Delta_C_bat = I_bat * t_Flug + Delta_C_bat;                 % entnommene Ladung
C_Rest_V = (C_bat_Peukert - Delta_C_bat) / C_bat_Peukert;   % Ladezustand als Verhaltnis

end