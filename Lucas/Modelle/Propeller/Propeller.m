function [Omega_neu,tau] = Propeller(V_A, alpha, Thrust, RPM_map, V_map, T_map, Tau_map)

% PROPELLER Drehzahl und Drehmomentberechnung

%       [omega_neu,tau] = Propeller(vi0, V_A, alpha, Thrust) ermittelt aus 
%       einem gegebenen Propellerkennfeld die Drehzahl für einen benötigten 
%       Schub sowie das sich daraus ergebende Drehmoment durch lineare 
%       Interpolation.



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

Omega_neu = rpm/60*2*pi;            


% DREHMOMENT (analog zur Drehzahl)

Tau_interp = zeros(length(RPM_map),1);
tau_row = 0;


% 1. Interpolation

for i = 1:length(RPM_map)                                 % Berechnung des interpolierten Wertes für x_val in jeder Zeile
    tau_row = Tau_map(i,:);
    Tau_interp(i) = interp1(V_map,tau_row,v);
end

% 2. Interpolation

tau = interp1(RPM_map,Tau_interp,rpm);


end

