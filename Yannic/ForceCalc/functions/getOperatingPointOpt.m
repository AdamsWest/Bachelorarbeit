function [aero,prop,motor,esc,bat,env,mission] = ...
    getOperatingPointOpt(aero,prop,motor,esc,bat,env,mission)


i_length = length(mission.params.V_K{mission.i_configuration});
j_length = length(mission.params.gamma{mission.i_configuration});
n = i_length * j_length;
aero_tmp = initAerodynamicsProperties(aero, n);
prop_tmp = initPropellerProperties(prop, n);
motor_tmp = initMotorProperties(motor, n);
esc_tmp = initEscProperties(esc, n);
bat_tmp = initBatteryProperties(bat, n);
env_tmp = initEnvironmentProperties(env, n);
mission_tmp = initMissionProperties(mission, n);

i_iter = 1;

for i = 1:i_length
    
    for j = 1:j_length       
        
        % set the parameters to the mission
        mission.V_K = mission.params.V_K{mission.i_configuration}(i);
        mission.gamma = mission.params.gamma{mission.i_configuration}(j);

        % get the operating point of the aircraft
        [aero,prop,motor,esc,bat,env,mission] = ...
            getOperatingPoint(aero,prop,motor,esc,bat,env,mission);

        % store the results in the _time class
        aero_tmp = storeAerodynamicsState( aero_tmp, aero, i_iter );
        prop_tmp = storePropellerState( prop_tmp, prop, i_iter );
        motor_tmp = storeMotorState( motor_tmp, motor, i_iter );
        esc_tmp = storeEscState( esc_tmp, esc, i_iter );
        bat_tmp = storeBatteryState( bat_tmp, bat, i_iter );
        env_tmp = storeEnvironmentState( env_tmp, env, i_iter );
        mission_tmp = storeMissionState( mission_tmp, mission, i_iter );

        i_iter = i_iter + 1;
        
    end

end

% check whether limits are reached
isLimitReached = ( max(prop_tmp.workload,[],2) > 1 | ...
    max(motor_tmp.workload,[],2) > 1 | ...
    max(bat_tmp.workload,[],2) > 1 );

% set C rate to NaN if limits are reached
bat_tmp.C_rate(isLimitReached) = NaN;
% compute the number of possible operating points
n_possible = sum( not( isLimitReached ) );
% if there is no possible operating point, break
if n_possible == 0
    disp('No operating point could be found');
    return;
end

% determine the most efficient operating point
switch mission.params.configuration{mission.i_configuration}
    case 'loiter'
        [~,i_iter_best] = min( bat_tmp.C_rate );
    case 'climb'
        [~,i_iter_best] = min( bat_tmp.C_rate ./ ...
            ( mission_tmp.V_K .* sin(mission_tmp.gamma) ) );
    case 'cruise'
        [~,i_iter_best] = min( bat_tmp.C_rate ./ ...
            ( mission_tmp.V_K .* cos(mission_tmp.gamma) ) );
end

% pick the optimum result/operating point
aero = pickAerodynamicsState( aero_tmp, i_iter_best );
prop = pickPropellerState( prop_tmp, i_iter_best );
motor = pickMotorState( motor_tmp, i_iter_best );
esc = pickEscState( esc_tmp, i_iter_best );
bat = pickBatteryState( bat_tmp, i_iter_best );
env = pickEnvironmentState( env_tmp, i_iter_best );
mission = pickMissionState( mission_tmp, i_iter_best );

    
end