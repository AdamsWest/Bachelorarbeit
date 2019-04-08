function [aero,prop,motor,esc,bat,env,mission] = createProperties( ...
    settingsFileName, DATA_APC,axi_motor_db,Elektromodellflug )

% load the defined parameters from the file
run( settingsFileName );


%% set the aerodynamic class and properties
% initialize the aerodynamics class
aero = AerodynamicsClass;
% set the parameters in the params struct of the class
aero.params = aeroParams;
% set the mass of the aircraft
aero.params.m = aeroParams.m_empty + propParams.n*motorParams.m + ...
    aeroParams.payload + batParams.m;


%% set the propeller class and properties
% initialize the propeller class
prop = PropellerClass;
% set the parameters in the params struct of the class
prop.params = propParams;
% compute the propeller map from the parameters and set it to the params
[pmap.RPM, pmap.V, pmap.T_0, pmap.P_0, pmap.tau_0] = ...
    propellerMap( DATA_APC, prop.params.name );
prop.params.map = pmap;
prop.params.Omega_max = max(prop.params.map.RPM) * 2*pi/60;
% initialize the propeller map depending on the air density
prop = setPropellerMapRho(prop,1,1);
% compute geometric parameters
prop = setPropGeometicParams(prop);


%% set the motor class and properties
% initialize the motor class
motor = MotorClass;
% set the parameters in the params struct of the class
motor.params = motorParams;
% load specific motor parameters
if ~motorParams.isCustom
    motor = setSpecificMotorParams(motor,'axi_motor_db');
end


%% set the ESC class and properties
% initialize the esc class
esc = EscClass;
% set the parameters in the params struct of the class
% esc.params = escParams;


%% set the battery class and properties
% initialize the battery class
bat = BatteryClass;
% set the parameters in the params struct of the class
bat.params = batParams;
% compute battery mass or capacity
bat = computeBatteryMassOrCapacity(bat);
% set parameters for dynamic voltage computation if desired
if bat.params.isDynVolt
    % get typical points on the discharge curve for LiPos
    [ p.SoC_full, p.SoC_nom, p.SoC_exp, p.V_full, p.V_exp, p.V_nom, ...
        p.C_rate ] = batteryAverageParams( Elektromodellflug );
    % get the resistance according to a model
    Q = p.SoC_full * bat.params.C/3600;
   	p.R = 0.1077 / (Q*(0.1555*Q+0.9825*bat.params.C_rate_max)^0.5485);
    p.R = 0.0158;
    p.R_over_C = p.R / ( bat.params.C / 3600 );
    bat.params.dynVolt = p;
    % get the parameters of the discharge curve
    [bat.params.dynVolt.Eo, bat.params.dynVolt.A, ...
        bat.params.dynVolt.K] = batteryDischargeParams( ...
        [p.SoC_full,p.SoC_nom,p.SoC_exp,p.V_full,p.V_exp,p.V_nom,p.C_rate,p.R/(bat.params.C/3600)] );
    bat.params.dynVolt.B = 3/p.SoC_exp;
end



%% set the environment class and properties
% initialize the environment class
env = EnvironmentClass;
% set the parameters in the params struct of the class
env.params = envParams;
% compute the parameters at the stratopause
env = stratopauseConditions(env);


%% set the mission class and properties
%
mission = MissionClass;
%type
mission.params = missionParams;
%
mission.params.n_missionElements = length(missionParams.configuration);
%
if mission.params.isDynamic
    mission.params.n_operatingPoints = floor( mission.params.maxTime / ...
        mission.params.dt );
else
    mission.params.n_operatingPoints = mission.params.n_missionElements;
end
%
mission.params.H_0 = env.params.H_0;



%% initialize properties of the classes

n = mission.params.n_operatingPoints + 1;
aero = initAerodynamicsProperties( aero, n );
prop = initPropellerProperties( prop, n );
motor = initMotorProperties( motor, n );
esc = initEscProperties( esc, n );
bat = initBatteryProperties( bat, n );
env = initEnvironmentProperties( env, n );
mission = initMissionProperties( mission, n );


end