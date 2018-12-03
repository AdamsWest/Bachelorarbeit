%% Batteriemodell
clear
close all
clc

load('Elektromodellflug.mat');

id_bat = 26;
PWM = 0.80;
eta_PWM = 10.7;
I_mot = 10;
n_Prop = 4;

BDD_b = evalin('base','Elektromodellflug');
[Eo, A, K] = Batterie_parameter(cell2mat(BDD_b(id_bat,3)));  

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
hold on


%% for comparison

Q = 1.0618;                   % Q
Qnom = 0.9102;                % Qnom
Qexp = 0.2083;                % Qexp
Vfull = 4.0399;               % Vfull
Vexp = 3.7783;                % Vexp
Vnom = 3.5181;                % Vnom
i = 29.7754;                  % i
R = 0.0142;               % R_bat
B = 3/Qexp;
Batterie_data = [Q Qnom Qexp Vfull Vexp Vnom i R 0 0 0];        % Zwischenbelegung
[Eo,A,K] = Batterie_parameter(Batterie_data);
Batterie_data = [Q Qnom Qexp Vfull Vexp Vnom i R Eo A K];       % vollständiger Vektor


I_bat = PWM * I_mot / eta_PWM * n_Prop;

i_int = 0;
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
plot(x,V_bat,'r-')



