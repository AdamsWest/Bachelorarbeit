function [aero,prop,motor,esc,bat,env,mission] = ...
    getOperatingPoint(aero,prop,motor,esc,bat,env,mission)
% getOperatingPoint computes the operating points of the whole aircraft
%
%   The operating points of all component, aerodynamic, propeller, motor,
%   ESC, battery, environment and the mission state, is computed in this
%   function.
%
% Syntax:   [aero,prop,motor,esc,bat,env,mission] = ...
%               getOperatingPoint(aero,prop,motor,esc,bat,env,mission);
%
% Inputs:
%   aero        object of class AerodynamicsClass
%   prop        object of class PropellerClass
%   motor       object of class MotorClass
%   esc         object of class EscClass
%   bat         object of class BatteryClass
%   env         object of class EnvironmentClass
%   mission     object of class MissionClass
%
% Outputs:
%   aero        object of class AerodynamicsClass
%   prop        object of class PropellerClass
%   motor       object of class MotorClass
%   esc         object of class EscClass
%   bat         object of class BatteryClass
%   env         object of class EnvironmentClass
%   mission     object of class MissionClass
%
%
% See also: getOperatingPointOpt
%
%   Copyright 2019 TU Braunschweig
% *************************************************************************

%% environment
%
h_mid = mission.V_K * sin(mission.gamma) / 2 + mission. altitude;
[env.rho,env.a,env.p,env.T] = atmosphere(h_mid,env.params);
env.g = env.params.g;
%% aerodynamics
% 
ap = aero.params;
% check whether the aircraft is a multicopter or an airplane
if aero.params.type == 0
    % multicopter
    [ aero.Thrust, aero.Theta_1, aero.V_A, aero.alpha ] = ...
        aerodynamicsMulticopter( mission.V_K, mission.gamma, ...
        ap.c_W_copter_sideways, ap.c_W_copter_upper, ap.c_A_copter_max, ...
        ap.A_copter, ap.m, env.g, env.u_Wg, env.rho );
    % compute the vertical velocity component
    V_A_Prop = - aero.V_A * sin( aero.alpha );
else
    % airplane
    [ aero.Thrust, aero.V_A, ~ ] = ...
        aerodynamicsAirplane( ap.m, env.g, ap.E_star, ...
        ap.V_star, ap.rho_star, mission.gamma, env.rho );
    % compute the vertical velocity component
    V_A_Prop = aero.V_A;
    mission.V_K = aero.V_A;
    aero.alpha = 0;
end

%% propeller     
%
prop.T = aero.Thrust;
% compute the propeller operating point
[ prop.Omega, prop.tau ] = propellerOp( V_A_Prop, prop.T, ...
    prop.params.map.RPM, prop.params.map.V, prop.mapRho.T, prop.mapRho.tau );
% compute the mach number of the blade tips
V_tip = prop.Omega * prop.params.r;
prop.M_tip = V_tip / env.a;
% compute propeller induced velocity and efficiency
[ vi, mu_z, ~ ] = propellerInducedVelocity( aero.V_A, aero.alpha, ...
    aero.Thrust, env.g, env.rho, prop.params.F, prop.params.n );
prop.eta = prop.T * ( mu_z + vi ) / ( prop.Omega * prop.tau );   
% compute the relation of physical quantaties and technical limit
prop = getPropellerWorkload(prop);

%% motor
% compute the motor operating point
[ motor.U, motor.I, motor.eta ] = motorOp( prop.tau, motor.params.K_V, ...
    motor.params.I_0, motor.params.R_i, prop.Omega);
% compute the relation of physical quantaties and technical limit
motor = getMotorWorkload(motor);

%% ESC
% compute the ESC state
[ esc.PWM, esc.eta, esc.I_bat ] = escOp( motor.U, bat.U, motor.I, ...
    prop.params.n );

%% battery
% compute the battery C-rate, 1/h
bat.C_rate = esc.I_bat / bat.params.C * 3600;
% compute the battery voltage, V
bat.U = batteryVoltage(bat.SoC,bat.C_rate,bat.params.dynVolt) * ...
    bat.params.n_cell;
% compute the battery workload, -
bat = getBatteryWorkload(bat);
            
end
            