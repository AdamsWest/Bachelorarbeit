function [Thrust,Theta_1,V_A,alpha] = Aerodynamik(u_Wg,V_Kg,gamma,c_W_copter_seitlich,c_W_copter_oben,c_A_max,rho,A_copter,m,g)

% AERODYNAMIK   Berechnet die Aerodynamik eines Multicopters
%  
%   Auf Grundlage eines einfachen aerodynamischen Modells wird mithilfe der
%   Bahngeschwindigkeit, den dimensionslosen Beiwerten und den Copter
%   Spezifikationen der benötigten Schub und die Fluggeschwindigkeit berechnet.
%
% Syntax:  [Thrust,Theta_1,V_A,alpha] = Aerodynamik(u_Wg,V_Kg,gamma,c_W_copter_seitlich,c_W_copter_oben,c_A_max,rho,A_copter,m,g)
%
% Inputs:
%   u_Wg     Seitenwindgeschwindigkeit
%   V_Kg     Bahngeschwindigkeit
%   gamma    Bahnneigungswinkel
%   c_W_copter_seitlich    seitlicher, dimensionsloser Widerstansbeiwert
%   c_W_copter_oben        oberer, dimensionsloser Widerstandsbeiwert
%   c_A_max  maximaler Auftriebsbeiwert des Quadrocopters (bei +/-45° Anstellwinkel)
%   rho      Luftdichte
%   A_copter obere Stirnflaeche 
%   m        Masse
%   g        Erdbeschleunigung
%
% Outputs:
%   Thrust   benötigter Schub
%   Theta_1  Neigungswinkel
%   V_A      Absolute Fluggeschwindigkeit
%   alpha    Anstellwinkel
%
% Example:
%    Line 1 of example
%    Line 2 of example
%    [ values, derivatives ] = Untitled( x, myStruct.y, 2 )
%    [ values, ~ ] = Untitled( [ 1:.1:100 ] , myStruct.y, 1 )
%
% See also: PROPELLER,  MOTOR

%   Copyright 2018 TU-Braunschweig
% ******************************************************************************  

% Wenn kein Seitenwind vorhanden, setze ihn auf 0.01
if u_Wg == 0        
    u_Wg = 0.01;
end

% Fluggeschwindigkeit in x_g- und z_g-Richtung
w_Kg = - V_Kg * sin(gamma);
u_Kg = V_Kg * cos(gamma);

V_A = sqrt((u_Kg + u_Wg)^2 + w_Kg^2);       % Absolute Flugwindgeschwindigkeit
gamma_a = atan(- w_Kg / (u_Kg + u_Wg));     % Windanstellwinkel

% Initialisierungen
Theta_1 = 0;
Delta_Theta = 2;
i = 0;
    while Delta_Theta > 0.001*pi/180                % Abbruchkriterium
        alpha = - gamma_a + Theta_1;                % Berechne alpha
        c_W = - (c_W_copter_oben - c_W_copter_seitlich)/2 * cos(2*alpha) + (c_W_copter_oben + c_W_copter_seitlich)/2; % Berechne den Widerstandsbeiwert in Abhängigkeit von alpha
        c_A = c_A_max * sin(2*alpha);               % Berechne den Auftriebsbeiwert in Abhängigkeit von alpha
        W = c_W * rho/2 * V_A^2 * A_copter;         % Widerstandskraft
        A = c_A * rho/2 * V_A^2 * A_copter;         % Auftriebskraft
        X_g = - W * cos(gamma_a) - A * sin(gamma_a);        % resultierende Kraftkomponente in x_g_Richtung
        Z_g = W * sin(gamma_a) - A * cos(gamma_a) + m*g;    % resultierende Kraftkomponente in z_g_Richtung
        Theta_2 = -atan(-X_g / Z_g);                % numerische Annäherung des Neigungswinkels
        Delta_Theta = abs(Theta_2 - Theta_1);
        Theta_1 = Theta_2;  
        i = i + 1;
    end
    Thrust = sqrt(X_g^2 + Z_g^2);                   % Berechne Schub aus den Kraftkomponenten mithilfe des Satz des Pythagoras
    
end