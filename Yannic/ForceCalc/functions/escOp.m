function [ PWM, eta_PWM, I_bat ] = escOp( U_mot, U_bat, I_mot, n_Prop )
% escOp computes the operating point of the motor controller (ESC)
%
%   The operating point of the ESC is computed based on the model of
%   Lubrano (2016). It considers an efficiency factor (eta) which depends
%   on the duty cycle of the ESC (PWM).
%
% Syntax:   [ PWM, eta_PWM, I_bat ] = escOp( U_mot, U_bat, I_mot, n_Prop )
%
% Inputs:
%   U_mot       motor voltage (scalar), in V
%   U_bat       battery voltage (scalar), in V
%   I_bat       battery current (scalar), in A
%
% Outputs:
%   PWM         pulse width modulation or duty cycle (scalar), -
%   eta_PWM     efficiency of the ESC (scalar), -
%
%
% See also: batteryVoltage
%
%   Copyright 2019 TU Braunschweig
% *************************************************************************

% computation of the duty cycle
PWM = U_mot / U_bat;
       
% computation of the efficiency
if PWM > 0 && PWM < 0.5 
    eta_PWM = 0.7 * PWM + 0.50;
elseif PWM >= 0.5 && PWM <= 1
	eta_PWM = 0.2 * PWM + 0.75;
else
	eta_PWM = 1;
end

% computation of the battery current
I_bat = PWM * I_mot / eta_PWM * n_Prop;

end