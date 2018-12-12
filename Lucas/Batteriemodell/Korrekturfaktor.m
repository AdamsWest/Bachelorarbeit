%% Batteriemodell
clear
close all
clc

%% Allgemeine Parameter
load('Elektromodellflug.mat');



% Anmerkung: id_bat 33 zu geringe Spannung
PWM = 0.80;
eta_PWM = 0.7;
I_mot = 8;
n_Prop = 6;
%% Initialisierung für Fehlertoleranz
control =0;
Toleranz = zeros(length(Elektromodellflug),1);


for n = 1:length(Elektromodellflug)
    
    id_bat = n;
    
    
    
    %% set time step delta t in seconds
    step_dt = 1;
    
    BDD_b = evalin('base','Elektromodellflug');
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
    V_bat_1 = zeros(2800,1);
    control = 0;
    bar = zeros(2800,1);
    
    I_bat = PWM * I_mot / eta_PWM * n_Prop;
    
    for i = 1:2800
        
        
        i_int = I_bat*step_dt + i_int;   % integral of the current
        
        
        % calculating the battery voltage (of one cell)
        V_bat_1(i) = Eo - R*I_bat - K * Q / (Q - i_int) * (i_int + I_bat*0) + A * exp(-B*i_int);
        V_bat_1(i) = N_el * V_bat_1(i);   % the battery voltage of all cells
        
        
%         if V_bat_1 < (3.1 * N_el)
%             V_bat_1 = NaN;
%             % break
%         end
        
        if V_bat_1(i) < 3.1*N_el
             control = 1;
        end
        if control == 1
            V_bat_1(i) = NaN;
        end
        
        bar(i) = 3.1*N_el;
    end
    % Kurve der Batterie plotten

    
    
    
    
    
    
    
    %% for comparison
    %%
    % insert capacity later if needed
    Cnom = Elektromodellflug{id_bat,5}/1000;                        % Cnom
    Q = 1.0618;                                                     % Q
    Qnom = 0.9102;                                                  % Qnom
    Qexp = 0.2083;                                                  % Qexp
    Vfull = 4.0399;                                                 % Vfull
    Vexp = 3.7783;                                                  % Vexp
    Vnom = 3.5181;                                                  % Vnom
    i = 1/100;                                                      % i
    R = 0.015;                                                      % R_bat
    B = 3/Qexp;
    Batterie_data = [Q Qnom Qexp Vfull Vexp Vnom i R 0 0 0];        % Zwischenbelegung
    [Eo,A,K] = Batterie_parameter(Batterie_data);
    %Batterie_data = [Q Qnom Qexp Vfull Vexp Vnom i R Eo A K];       % vollständiger Vektor
    
    
    I_bat = PWM * I_mot / eta_PWM * n_Prop;                         % Batteriestrom
    I_bat = I_bat/(Cnom);                                           % Normierung des Batteriestroms
    i_int = 0;
    V_bat_2 = zeros(2800,1);
    control = 0;
    
    for i = 1:2800
        
        
        i_int = I_bat*step_dt/3600 + i_int;   % integral of the current
        
        
        % calculating the battery voltage (of one cell)
        V_bat_2(i) = Eo - R*I_bat - K * Q / (Q - i_int) * (i_int + I_bat*0) + A * exp(-B*i_int);
        % the battery voltage of all cells
        V_bat_2(i) = N_el * V_bat_2(i)*1.065;  % <-- KORREKTURFAKTOR VON DURCHSCHNITTLICH 5%
        
        
        if V_bat_2(i) < 3.1*N_el
             control = 1;
        end
        if control == 1
            V_bat_2(i) = NaN;
        end
        
        
    end
    
    for i = 2800:-1:1
        if isnan(V_bat_1(i)) == 1
            V_bat_1(i) = 1;
        end
        if isnan(V_bat_2(i)) == 1
            V_bat_2(i) = 1;
        end
    end
%     V_bat_2 = V_bat_2(1:length(V_bat_1));
%     if isempty(V_bat_1) == 1 && isempty(V_bat_2) == 1
%         V_bat_1 = 1;
%         V_bat_2 = 1;
%     end
    
    
    Toleranz(n) = mean(abs((V_bat_2)./V_bat_1-1)*100);
    
end




%% Fehlertoleranz
x = 1:length(Elektromodellflug);
plot(x,Toleranz)
xlabel('Batterienummer');
ylabel('Durchschnittliche Abweichung in %');
Durchschnitt_Toleranz = mean(Toleranz);
