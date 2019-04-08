function U_cell = batteryVoltage(SoC,C_rate,params)
% batteryVoltage computes the voltage of the battery
%   
%   The battery voltage  of one cell is computed according to a discharge 
%   curve depending on the state of charge (SoC) and the discharge rate 
%   (C rate) of the battery.
%
% Syntax:  U_cell = batteryVoltage(SoC,C_rate,params)
%
% Inputs:
%   SoC         state of charge (scalar), -
%   C_rate      discharge rate (scalar), in 1/h
%   params      a struct containing several parameters
%
% Outputs:
%   U_cell      the voltage of one cell (scalar), in V
%
%
% See also: batteryDischargeParams, batteryAverageParams
%
%   Copyright 2019 TU Braunschweig
% *************************************************************************

% get parameters
Eo = params.Eo;
SoC_full = params.SoC_full;
R_over_C = params.R_over_C;
K = params.K;
A = params.A;
B = params.B;
% compute voltage
U_cell = Eo - R_over_C*C_rate - K*SoC_full ./ ( SoC_full -  (1- SoC ) ) .* ...
    ((1-SoC) + C_rate*0) + A * exp(-B*(1-SoC));

end