function [Omega,tau] = Propeller(V_A, alpha, Thrust, RPM_map, V_map, T_map, TAU_map)

% PROPELLER Interpolation der Drehzahl und des Drehmoments
%   
%   Mittels linearer Interpolation wird für einen benötigten Schub und die Fluggeschwindigkeit
%   aus einem gegebenen Propellerkennfeld die Drehzahl und das sich daraus ergebende Drehmoment 
%   ermittelt.
%
% Syntax:  [Omega,tau] = Propeller(V_A, alpha, Thrust, RPM_map, V_map, T_map, TAU_map)
%
% Inputs:
%    V_A	absolute Fluggeschwindigkeit
%    alpha	Rotoranstellwinkel
%    Thrust	benötigter Schub
%    RPM_map	Drehzahlkennfeld
%    V_map	Geschwindigkeitskennfeld
%    T_map	Schubkennfeld
%    TAU_ma	Drehmomentkennfeld
%
% Outputs:
%    Omega	Drehzahl
%    tau	Drehmoment
%
% Example:
%    Line 1 of example
%    Line 2 of example
%    [ values, derivatives ] = Untitled( x, myStruct.y, 2 )
%    [ values, ~ ] = Untitled( [ 1:.1:100 ] , myStruct.y, 1 )
%
% See also: AERODYNAMIK,  MOTOR

%   Copyright 2018 TU-Braunschweig
% ******************************************************************************


% Berechnung der Drehzahl und des Drehmoments aus dem Propellerkennfeld

v = - V_A * sin(alpha);                                   % Berechnung der Durchflussgeschwindigkeit senkrecht durch die Rotorebene


% DREHZAHL


T_interp = zeros(length(RPM_map),1);                      % Initialisierung eines Ergebnisvektors der ersten interpolierten Werte
t_row = 0;                                                % Initialisierung der Spalten


% 1. Interpolation

for i = 1:length(RPM_map)                                 % Berechnung des interpolierten Wertes für x_val in jeder Zeile
    t_row = T_map(i,:);
    T_interp(i) = interp1(V_map,t_row,v);
end

% 2. Interpolation

rpm = interp1(T_interp,RPM_map,Thrust);                   % Bestimmung der Drehzahl in (U/min)

Omega = rpm/60*2*pi;    				  % Angabe der Drehzahl in U/s	        


% DREHMOMENT (analog zur Drehzahl)

Tau_interp = zeros(length(RPM_map),1);			  % Initialisierung eines Ergebnisvektors der ersten interpolierten Werte
tau_row = 0;						  % Initialisierung der Spalten


% 1. Interpolation

for i = 1:length(RPM_map)                                 % Berechnung des interpolierten Wertes für x_val in jeder Zeile
    tau_row = TAU_map(i,:);
    Tau_interp(i) = interp1(V_map,tau_row,v);
end

% 2. Interpolation

tau = interp1(RPM_map,Tau_interp,rpm);			  % Interpolation des endgültigen Drehmoments in Nm


end

