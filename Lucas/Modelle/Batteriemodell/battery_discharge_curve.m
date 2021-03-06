%function time = battery_discharge_curve(BDD_b,id_bat,i_m,v_m,V_bat_min,Eo,A,K,i_int,step_dt)

function [V_bat] = battery_discharge_curve(i_m,v_m,Eo,A,K,Q,Qnom,Qexp,Vfull,Vnom,Vexp,R,C_Bat,N_el,i_int,step_dt,i_bat)

% battery parameters
%Q = BDD_b{id_bat,3}(1);
B = 3/Qexp;
%R = BDD_b{id_bat,3}(8);
i_star = 0;
%N_el = BDD_b{id_bat,4};
%C = BDD_b{id_bat,5} / 1000;
% C_rate_max = BDD_b{id_bat,6}; wird au�erhalb Funktion ben�tigt


% initializing the variables
% i_int = zeros(3,1));
i_int = i_bat*step_dt + i_int;   % integral of the current
% calculating the battery voltage (of one cell)
V_bat = Eo - R*i_bat - K * Q / (Q - i_int) * ...
    (i_int + i_bat*i_star) + A * exp(-B*i_int);
V_bat = N_el * V_bat;   % the battery voltage of all cells
end
