function temps_vol = batterie_curbe_de_decharge(BDD_b,id_bat,i_m,v_m,V_bat_min,Eo,A,K)


%% set time step delta t in seconds
step_dt = 1; 
%%

% battery parameters
Q = BDD_b{id_bat,3}(1);
B = 3/BDD_b{id_bat,3}(3);
R = BDD_b{id_bat,3}(8);
N_el = BDD_b{id_bat,4};
C = BDD_b{id_bat,5} / 1000;
C_rate_max = BDD_b{id_bat,6};

% initializing the variables
i_int = zeros(3,length(i_m));
V_bat = zeros(3,length(i_m));
T = ones(1,length(i_m))*2;
temps_vol = zeros(1,length(i_m));

%calculating the battery state
for k = 1:length(i_m)           % for every combination of motors and propellers
    
    % initializating the variables
    eta_pwm = ones(2,1);
    PWM = ones(2,1);
    i_bat = ones(2,1)*i_m(k);
    
    % using the formular of Trembley
    while isnan(V_bat(T(k)-1,k)) ~= 1   % while V_bat is not NaN 
        i_int(T(k)+1,k) = i_bat(T(k),1)*step_dt + i_int(T(k),k);   % integral of the current 
        % calculating the battery voltage (of one cell)
        V_bat(T(k),k) = Eo - R*i_bat(T(k),1) - K * Q / (Q - i_int(T(k)+1,k)) * (i_int(T(k)+1,k) + i_bat(T(k),1)*0) + A * exp(-B*i_int(T(k)+1,k));
        V_bat(T(k),k) = N_el * V_bat(T(k),k);   % the battery voltage of all cells
        
        % if the battery voltage is lower than the minimum allowable
        % voltage, set the voltage to NaN (to exit the calculation in the
        % following call of the while loop
        if V_bat(T(k),k) < V_bat_min * N_el
            V_bat(T(k):end,k) = NaN;
            break
        end  
        
        % calculating the PWM, the PWM efficiency and the PWM loss
        PWM(T(k),1) = v_m(k) / V_bat(T(k),k);
        if PWM(T(k),1) <= 0.5
            eta_pwm(T(k),1) = 0.7 * PWM(T(k),1) + 0.5;
        elseif PWM(T(k),1) > 1
            eta_pwm(T(k),1) = 1;
        else
            eta_pwm(T(k),1) = 0.2 * PWM(T(k),1) + 0.75;
        end
        % this is the battery current considering the PWM loss
        i_bat(T(k)+1,1) = i_m(k) * v_m(k) / V_bat(T(k),k) / eta_pwm(T(k),1);
        
        % counting the time (in seconds)
        T(k) = T(k) + 1;
        if T > 86400   % break after a battery run-time of one day to not get stuck in the loop (in case of an error)
            break
        end
    end

    % if the discharge rate is higher than the maximum allowable discharge
    % rate, stop the calculation or set the flight time to zero
    if i_bat < C_rate_max * C
        temps_vol(k) = T(k) - 3;
    else
        temps_vol(k) = 0;
    end

end

temps_vol = temps_vol';

end
