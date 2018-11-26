%% Batteriemodell
clear
close all
clc



id_bat = 4;

BDD_b = evalin('base','Elektromodellflug');
[Eo, A, K] = batterie_parametres(cell2mat(BDD_b(id_bat,3)));  

%temps_vol = batterie_curbe_de_decharge(BDD_b,id_bat,i_m,v_m,V_bat_min,Eo,A,K)


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
T = 2;
temps_vol = 0;

I_bat = PWM * I_mot / eta_PWM * n_Prop;

while isnan(V_bat0) ~= 1   % while V_bat is not NaN
    
    
    i_int = I_bat*step_dt + i_int;   % integral of the current
    
    
    % calculating the battery voltage (of one cell)
    V_bat = Eo - R*I_bat - K * Q / (Q - i_int) * (i_int + I_bat*0) + A * exp(-B*i_int);
    V_bat = N_el * V_bat;   % the battery voltage of all cells
    
    T = T + 1;
    if T > 86400   % break after a battery run-time of one day to not get stuck in the loop (in case of an error)
        break
    end
end