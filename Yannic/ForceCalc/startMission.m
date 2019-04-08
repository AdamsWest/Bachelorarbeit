% initialization of all objects
initAircraft;

for i_time = 1:length(mission.time)-1
    
    % create instances with scalar properties of the objects
    aero_t = pickAerodynamicsState( aero, i_time );
    prop_t = pickPropellerState( prop, i_time );
    motor_t = pickMotorState( motor, i_time );
    esc_t = pickEscState( esc, i_time );
    bat_t = pickBatteryState( bat, i_time );
    env_t = pickEnvironmentState( env, i_time );
    mission_t = pickMissionState( mission, i_time ); 
    
    % get the optimum operating point
    [aero_t,prop_t,motor_t,esc_t,bat_t,env_t,mission_t] = ...
        getOperatingPointOpt(aero_t,prop_t,motor_t,esc_t,bat_t,env_t,mission_t);
    
%     % if limit is reached, set the remaining capacity to NaN
%     bat_t.C_rest( isLimitReached ) = NaN;
    
    % store the most efficient operating point
    aero = storeAerodynamicsState( aero, aero_t, i_time );
    prop = storePropellerState( prop, prop_t, i_time );
    motor = storeMotorState( motor, motor_t, i_time );
    esc = storeEscState( esc, esc_t, i_time );
    bat = storeBatteryState( bat, bat_t, i_time );
    env = storeEnvironmentState( env, env_t, i_time );
    mission = storeMissionState( mission, mission_t, i_time );
    
    % get flight time
    dt = getFlightTime( mission );
    
    % increment classes
    aero = incrementAerodynamics( aero, dt );
    prop = incrementPropeller( prop, dt );
    motor = incrementMotor( motor, dt );
    esc = incrementEsc( esc, dt );
    bat = incrementBattery( bat, dt );
    env = incrementEnvironment( env, dt );
    mission = incrementMission( mission, dt );
    
    % integrate the capacity from the current and the time step
    bat = integrateCapacity( bat );
    if bat.SoC(i_time+1) <= 0
        break;
    end
    
    % integrate the position from the velocity and the flight path angle
    mission = integratePosition( mission );
    
    % increment mission if mission part is completed
    mission = incrementConfiguration( mission );
    
    % stop computation if mission is completed
    if mission_t.i_configuration == 0
        disp('mission completed')
        break;
    end
    
    disp(['t=',num2str(mission.time(mission.idx)),'s'])
    
end

%% plot results

i = mission.idx-1;
%
subplot(4,3,1)
plot(mission.altitude(1:i),bat.SoC(1:i));
grid on
ylim([0, 1])
xlabel('altitude, m')
ylabel('state of charge, -')
%
hold on
subplot(4,3,2)
yyaxis left
plot(mission.altitude(1:i),prop.Omega(1:i));
hold on
plot(mission.altitude(1:i),max(prop.params.Omega_max*ones(1,i)),'k--')
ylim([0, Inf])
ylabel('propeller angular velocity, 1/s')
yyaxis right
plot(mission.altitude(1:i),prop.M_tip(1:i))
plot(mission.altitude(1:i),ones(1,i))
grid on
ylim([0, Inf])
xlabel('altitude, m')
ylabel('blade tip Mach number, -')
%
subplot(4,3,3)
yyaxis left
stairs(mission.altitude(1:i),prop.T(1:i));
grid on
ylim([0 max(prop.T)])
xlabel('altitude, m')
ylabel('propeller thrust, N')
yyaxis right
stairs(mission.altitude(1:i),prop.tau(1:i));
ylim([0 max(prop.tau)])
xlabel('altitude, m')
ylabel('propeller torque, Nm')
%
subplot(4,3,4)
stairs(mission.altitude(1:i),motor.I(1:i));
grid on
ylim([0 motor.params.I_max])
xlabel('altitude, m')
ylabel('motor current, A')
%
subplot(4,3,5)
stairs(mission.altitude(1:i),motor.U(1:i));
grid on
ylim([0 max(bat.U)])
xlabel('altitude, m')
ylabel('motor voltage, V')
%
subplot(4,3,6)
stairs(mission.altitude(1:i),bat.C_rate(1:i));
grid on
ylim([0 bat.params.C_rate_max])
xlabel('altitude, m')
ylabel('battery C rate, 1/h')
%
subplot(4,3,7)
stairs(mission.altitude(1:i),bat.U(1:i));
grid on
ylim([0 max(bat.U)])
xlabel('altitude, m')
ylabel('battery voltage, V')
%
subplot(4,3,8)
stairs(mission.altitude(1:i),esc.PWM(1:i));
grid on
ylim([0 1])
xlabel('altitude, m')
ylabel('ESC duty cycle, -')
%
subplot(4,3,9)
stairs(mission.altitude(1:i),prop.eta(1:i).*motor.eta(1:i).*esc.eta(1:i));
hold on
stairs(mission.altitude(1:i),prop.eta(1:i));
stairs(mission.altitude(1:i),motor.eta(1:i));
stairs(mission.altitude(1:i),esc.eta(1:i));
grid on
ylim([0 1])
xlabel('altitude, m')
ylabel('efficiency, -')
legend('total','propeller','motor','ESC')
%
subplot(4,3,10)
stairs(mission.altitude(1:i),mission.gamma(1:i)*180/pi);
grid on
ylim([0 90])
xlabel('altitude, m')
ylabel('flight path angle, °')
%
subplot(4,3,11)
yyaxis left
stairs(mission.altitude(1:i),mission.V_K(1:i));
grid on
ylim([0 Inf])
xlabel('altitude, m')
ylabel('flight path velocity, m/s')
yyaxis right
stairs(mission.altitude(1:i),mission.V_K(1:i).*sin(mission.gamma(1:i)));
ylim([0 Inf])
xlabel('altitude, m')
ylabel('climb rate, m/s')
%
subplot(4,3,12)
stairs(mission.altitude(1:i),mission.time(1:i));
grid on
ylim([0 Inf])
xlabel('altitude, m')
ylabel('flight time, s')