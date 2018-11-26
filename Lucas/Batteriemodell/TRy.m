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
i_int = 0;
V_bat = 0;
T = 1;
temps_vol = 0;

I_bat = PWM * i_m / eta_PWM * n_Prop;

V_bat = zeros(1,100000);

while isnan(V_bat) ~= 1   % while V_bat is not NaN
    
    
    i_int = I_bat*step_dt + i_int;   % integral of the current
    
    
    % calculating the battery voltage (of one cell)
    V_bat(T) = Eo - R*I_bat - K * Q / (Q - i_int) * (i_int + I_bat*0) + A * exp(-B*i_int);
    V_bat(T) = N_el * V_bat(T);   % the battery voltage of all cells
    
    
    T = T + 1;
    if T > 86400   % break after a battery run-time of one day to not get stuck in the loop (in case of an error)
        break
    end
end
