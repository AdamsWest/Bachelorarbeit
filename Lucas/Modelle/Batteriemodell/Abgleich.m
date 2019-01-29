%% Batteriemodell
clear
close all
clc

load('Elektromodellflug.mat');


PWM = 0.80;
eta_PWM = 0.7;
I_mot = 10;
n_Prop = 4;

BDD_b = evalin('base','Elektromodellflug');
 



%% set time step delta t in seconds
step_dt = 1; 



%for m = 1:length(Elektromodellflug)
    id_bat = 26;
    [Eo, A, K] = Batterie_parameter(cell2mat(BDD_b(id_bat,3)));
    
    % battery parameters
    Q = BDD_b{id_bat,3}(1);
    B = 3/BDD_b{id_bat,3}(3);
    R = BDD_b{id_bat,3}(8);
    N_el = BDD_b{id_bat,4};
    C = BDD_b{id_bat,5} / 1000;
    C_rate_max = BDD_b{id_bat,6};
    
    % initializing the variables
    i_int = 0;
    V_bat = zeros(200,1);
    T = 2;
    temps_vol = 0;
    
    I_bat = PWM * I_mot / eta_PWM * n_Prop;
    
    for i = 1:2800
        
        
        i_int = I_bat*step_dt + i_int;   % integral of the current
        
        
        % calculating the battery voltage (of one cell)
        V_bat(i) = Eo - R*I_bat - K * Q / (Q - i_int) * (i_int + I_bat*0) + A * exp(-B*i_int);
        V_bat(i) = N_el * V_bat(i);   % the battery voltage of all cells
        
        
        if V_bat < 3.1 * N_el
            V_bat = NaN;
            break
        end
        
    end
    x = 1:2800;
    plot(x,V_bat)
    ylim([0 30]) 
    hold on
%end
%%
%% for comparison
%%

Cnom = 4;                     % Cnom
Q = 1.0618;                   % Q
Qnom = 0.9102;                % Qnom
Qexp = 0.2083;                % Qexp
Vfull = 4.0399;               % Vfull
Vexp = 3.7783;                % Vexp
Vnom = 3.5181;                % Vnom
i = 1/100;                    % i
R = 0.015;                   % R_bat
B = 3/Qexp;
Batterie_data = [Q Qnom Qexp Vfull Vexp Vnom i R 0 0 0];        % Zwischenbelegung
[Eo,A,K] = Batterie_parameter_2(Batterie_data);
Batterie_data = [Q Qnom Qexp Vfull Vexp Vnom i R Eo A K];       % vollständiger Vektor


I_bat = PWM * I_mot / eta_PWM * n_Prop;
I_bat = I_bat/(Cnom);

i_int = 0;
for i = 1:2800
    
    
    i_int = I_bat*step_dt/3600 + i_int;   % integral of the current
    
    
    % calculating the battery voltage (of one cell)
    V_bat(i) = Eo - R*I_bat - K * Q / (Q - i_int) * (i_int + I_bat*0) + A * exp(-B*i_int);
    V_bat(i) = N_el * V_bat(i);   % the battery voltage of all cells
    
    
    if V_bat < 3.1 * 6
        V_bat = NaN;
        break
    end
    
end
x = 1:2800;
plot(x,V_bat,'rx')



