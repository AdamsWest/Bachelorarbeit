function [U_mot,I_mot,Omega,alpha_75] = Antrieb(Thrust,Omega,vi0,V_A,alpha,R,Theta_75,rho,c_d0,a,I_mot,I_0,U_mot,K_V,R_i)

c_T_0 = Thrust / ( pi * R^4 * rho * Omega^2 );                  % Schubbeiwert statisch             
Phi_75_0 = atan( vi0 / ( 0.75*R*Omega ) );                      % induzierter Winkel statisch
c_Q_0 = rho/4 * (c_d0/2 + a*Phi_75_0*(Theta_75-Phi_75_0));      % Momentenbeiwert statisch
alpha_75_0 = (Theta_75 - Phi_75_0) * 180/pi;                    % Anstellwinkel bei 75% des Propellerradius statisch
azcusw = c_T_0 / (alpha_75_0*pi/180);


% Berechnung der induzierten Geschwindigkeit

v = vi0;
mu_z = -V_A*sin(alpha);
mu = V_A*cos(alpha);
krit = 1;
while krit > 0.0005
    f = v - mu_z - vi0^2 / sqrt(mu^2 + v^2);
    fs = 1 + v * vi0^2 / (mu^2 + v^2)^(3/2);
    v_i_neu = v - f/fs;
    krit = abs(v_i_neu - v) / v_i_neu;
    v = v_i_neu;
end
vi_vi0 = (v - mu_z) / vi0;


vi = vi0 * vi_vi0;                                              % induzierte Geschwindigkeit im stationaeren Steigflug


% Neuberechnung der Drehzahl

fun = @(Omega_find) azcusw * pi * R^4 * rho * ( Theta_75 - atan((vi+mu_z)/(Omega_find*R*0.75))) * Omega_find^2 - Thrust;
x0 = Omega*3;
options = optimoptions('fsolve','Display','off');
Omega_neu = fsolve(fun,x0,options);                             % Iteration gleichbedeutend mit while-Schleife


% Neuberechnung der Motor-Zustandsgroessen

Phi_75 = atan( (vi + mu_z) / ( 0.75*R*Omega_neu ) );            % induzierter Winkel
alpha_75 = (Theta_75 - Phi_75) * 180/pi;                        % Anstellwinkel bei 75% des Propellerradius
c_Q = rho/4 * (c_d0/2 + a*Phi_75*(Theta_75-Phi_75));            % Momentenbeiwert dynamisch
I_mot_neu = (c_Q * Omega_neu^2) / (c_Q_0 * Omega^2) * (I_mot - I_0) + I_0;      % resultierender Motorstrom    
%k = pi * rho * R^5;
%I_mot_neu = I_mot * (k * c_Q * Omega_neu^2 * K_V + I_0) / (k * c_Q_0 * Omega^2 * K_V + I_0);
U_mot = U_mot * (Omega_neu / K_V + R_i * I_mot_neu) / (Omega / K_V + R_i * I_mot);  % resultierende Motorspannung
I_mot = I_mot_neu;                                              % Wertuebergabe
Omega = Omega_neu;                                              % Wertuebergabe     


end






