function [PWM,eta_PWM] = ESC(U_mot,U_bat)

% ESC   Berechnet die Pulsweitenmodulation und den Wirkungsgrad des
%       Motorreglers
%
%   Mithilfe des Motorstroms und der nominellen Batteriespannung wird die
%   benötigte Pulsweitenmodulation und der Wirkungsgrad dieser nach dem
%   Modell von Lubrano(2016) berechnet.
%
% Syntax:  [PWM,eta_PWM] = ESC(U_mot,U_bat_nom)
%
% Inputs:
%   U_mot   Motorstrom
%   U_bat_nom   nominelle Batteriespannung
%
% Outputs:
%   PWM     Pulsweitenmodulation
%   eta_PWM Wirkungsgrad des Motorreglers
%
% Example:
%    Line 1 of example
%    Line 2 of example
%    [ values, derivatives ] = Untitled( x, myStruct.y, 2 )
%    [ values, ~ ] = Untitled( [ 1:.1:100 ] , myStruct.y, 1 )
%
% See also: BATTERIE, MOTOR

%   Copyright 2018 TU-Braunschweig
% ******************************************************************************

PWM = U_mot / U_bat;            % Pulsweitenmodulation
       
% Bestimmung des ESC-Wirkungsgrades

if PWM > 0 && PWM < 0.5 
    eta_PWM = 0.7 * PWM + 0.50;
elseif PWM >= 0.5 && PWM <= 1
	eta_PWM = 0.2 * PWM + 0.75;
else
	eta_PWM = 1;
end

eta_PWM = 1;
% eta_PWM = eta_PWM + (1-eta_PWM)/2;

end