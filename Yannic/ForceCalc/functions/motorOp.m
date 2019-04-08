function [ U_mot, I_mot, eta_mot ] = motorOp( tau, K_V, I_0, R_i, Omega )  
% motorOp computes the motor state
%
%   In a simple motor model [1] the motor state encompassing the motor 
%   voltage, the motor current and the motor efficiency is computated.
%
% Literature:
%   [1] Mark Drela. “First-Order DC Electric Motor Model”. Technical report
%       MIT, USA, 2007.
%
% Syntax:  [ U_mot, I_mot, eta_mot ] = motorOp( tau, K_V, I_0, R_i, Omega )
%
% Inputs:
%    tau        propeller torque (scalar), in Nm
%    K_V        motor Kv value (scalar), in 1/(V*s)
%    I_0        no load current (scalar), in A
%    R_i        internal motor resistance (scalar), in Ohm
%    Omega      propeller speed of rotation (scalar), in 1/s
%
% Outputs:
%    U_mot      motor voltage (scalar), in V
%    I_mot      motor current (scalar), in A
%    eta_mot    motor efficiency (scalar), in %
%
% See also: propellerOp,  escOp,  batteryVoltage

%   Copyright 2019 TU Braunschweig
% *************************************************************************

% Motor current
I_mot = tau*K_V + I_0;
% Motor voltage
U_mot = Omega/(K_V) + R_i*I_mot;            
% Motor efficiency
eta_mot = (tau*Omega)/(U_mot*I_mot);      

end

