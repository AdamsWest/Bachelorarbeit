function [v] = Propeller_umgk(Omega,PWM,V_map,RPM_map,T_map)

% PROPELLER berechnet die Fluggeschwindigkeit

%       [v] = Propeller_umgk(Omega,PWM,V_map,RPM_map,T_map) interpoliert die 
%       Fluggeschwindigkeit anhand der Drehzahl und des Schubes aus dem
%       Propellerkennfeld


rpm = Omega*60/(2*pi);

Thrust = PWM * max(max(T_map));

% 1. Interpolation

T_interp = zeros(1,length(V_map));
t_column = 0;

for i = 1:length(V_map)
    t_column = T_map(:,i);
    T_interp(i) = interp1(RPM_map,t_column,rpm);
end

% 2. Interpolation

v = interp1(T_interp,V_map,Thrust);

end

