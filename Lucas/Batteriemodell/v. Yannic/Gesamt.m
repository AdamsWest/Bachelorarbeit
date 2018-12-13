%% Untersuchung der Batterien %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear 
close all 
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% NORMZELLE ERZEUGEN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('Elektromodellflug.mat');
DATA = Elektromodellflug;

% Löschen der Ausreißer
DATA(14,:) = [];       % id_bat = 14
DATA(29,:) = [];       % id_bat = 30
DATA(38,:) = [];       % id_bat = 40
DATA(60,:) = [];       % id_bat = 63

% Initialisierungen
sum = 0;
sum_1 = 0;
sum_2 = 0;
sum_3 = 0;
sum_4 = 0;
sum_5 = 0;
sum_6 = 0;
% sum_7 = 0;
sum_8 = 0;


capacity = zeros(length(DATA),1);
resistance = zeros(length(DATA),1);

for i = 1:length(DATA)
    
    % Normierung von Q, Q_nom und Q_exp mit der Kapazität in As
    
    DATA{i,3}(1) = DATA{i,3}(1) * 1000 / (DATA{i,5}*3600);
    DATA{i,3}(2) = DATA{i,3}(2) * 1000 / (DATA{i,5}*3600);
    DATA{i,3}(3) = DATA{i,3}(3) * 1000 / (DATA{i,5}*3600);

    % Normzelle: (arithmetischer Mittelwert)
    
    sum_1 = sum_1 + DATA{i,3}(1);
    sum_2 = sum_2 + DATA{i,3}(2);
    sum_3 = sum_3 + DATA{i,3}(3);
    sum_4 = sum_4 + DATA{i,3}(4);
    sum_5 = sum_5 + DATA{i,3}(5);
    sum_6 = sum_6 + DATA{i,3}(6);
    % sum_7 = sum_7 + Elektromodellflug_norm{i,3}(7);
    sum_8 = sum_8 + DATA{i,3}(8);
    
    capacity(i) = DATA{i,5}/1000;
    resistance(i) = DATA{i,3}(8);
    
    % sum = sum + Elektromodellflug_norm{i,5}/(Elektromodellflug_norm{i,4}*3.6);

end


% arithmetischer Mittelwert über alle Batterien
% durchschnittliche Kapazität pro Zelle:
% Cnom = sum / length(Elektromodellflug_norm);
Q = sum_1 / length(DATA);
Qnom = sum_2 / length(DATA);
Qexp = sum_3 / length(DATA);
Vfull = sum_4 / length(DATA);
Vexp = sum_5 / length(DATA);
Vnom = sum_6 / length(DATA);
i = 1/100;
R = sum_8 / length(DATA);


% Batterieparameter
M_A = [1, 1, 0 ; 1, exp(-3), -Q/(Q-Qexp)*(Qexp+i) ; 1, exp(-3*Qnom/Qexp), -Q/(Q-Qnom)*(Qnom + i)];
b = [Vfull + R*i ; Vexp + R*i ; Vnom + R*i];
x = M_A\b;

Eo = x(1);
A = x(2);
K = x(3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FEHLERUNTERSUCHUNG %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



PWM = 0.80;
eta_PWM = 0.7;
I_mot = 8;
n_Prop = 6;
C_Rate = 15;

%% Initialisierung für Fehlertoleranz
control =0;
tolerance = zeros(length(DATA),1);
C_Rate_max = 50;                                        % Maximal zu untersuchende C_Rate festlegen
tolerance_crate = [1:length(DATA)]';       % Sammeln aller Abweichungen aller Batterien für jede C_Rate in 1er Schritten

for k = 1:1:C_Rate_max
    
    C_Rate = k;
    for n = 1:length(DATA)
        
        id_bat = n;
        
        
        
        %% set time step delta t in seconds
        step_dt = 1;
        

        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        %%  ORIGINALZELLE
        
        % DATA = evalin('base','Elektromodellflug');
        [Eo, A, K] = Batterie_parameter(cell2mat(DATA(id_bat,3)));
        % battery parameters
        Q = DATA{id_bat,3}(1);
        B = 3/DATA{id_bat,3}(3);
        R = DATA{id_bat,3}(8);
        N_el = DATA{id_bat,4};
        C = DATA{id_bat,5} / 1000;
        C_rate_max = DATA{id_bat,6};
        
        % initializing the variables
        i_int = 0;
        V_bat_1 = zeros(floor(3600*1.5/C_Rate),1);
        control = 0;
        bar = zeros(floor(3600*1.5/C_Rate),1);
        
        %     I_bat = PWM * I_mot / eta_PWM * n_Prop;
        I_bat = C_Rate*C;
        
        for i = 1:floor(3600*1.5/C_Rate)
            
            
            i_int = I_bat*step_dt + i_int;   % integral of the current
            
            
            % calculating the battery voltage (of one cell)
            V_bat_1(i) = Eo - R*I_bat - K * Q / (Q - i_int) * (i_int + I_bat*0) + A * exp(-B*i_int);
            V_bat_1(i) = N_el * V_bat_1(i);   % the battery voltage of all cells
            
     
            
            if V_bat_1(i) < 3.1*N_el
                control = 1;
            end
            if control == 1
                V_bat_1(i) = NaN;
            end
            
            bar(i) = 3.1*N_el;
        end
             
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        %% NORMZELLE
        Cnom = DATA{id_bat,5}/1000;
        i = 1/100;
        B = 3/Qexp;
        
        
        %% Vorgaben
        %     I_bat = PWM * I_mot / eta_PWM * n_Prop;                   % Batteriestrom
        I_bat = C*C_Rate;
        I_bat = I_bat/(Cnom);                                           % Normierung des Batteriestroms
        i_int = 0;
        V_bat_2 = zeros(floor(3600*1.5/C_Rate),1);
        control = 0;
        
        for i = 1:floor(3600*1.5/C_Rate)
            
            
            i_int = I_bat*step_dt/3600 + i_int;   % integral of the current
            
            
            % calculating the battery voltage (of one cell)
            V_bat_2(i) = Eo - R*I_bat - K * Q / (Q - i_int) * (i_int + I_bat*0) + A * exp(-B*i_int);
            % the battery voltage of all cells
            V_bat_2(i) = N_el * V_bat_2(i);  
            
            
            if V_bat_2(i) < 3.1*N_el
                control = 1;
            end
            if control == 1
                V_bat_2(i) = NaN;
            end
            
            
        end
        
        % Berechnung der Toleranz nur bis zum ersten Mal, wenn V_bat unter
        % V_bat < 3.1 * N_el, ansonsten alle Ergebnisse entfernen
        flaeche1 = trapz(V_bat_1(~isnan(V_bat_1)));
        flaeche2 = trapz(V_bat_2(~isnan(V_bat_2)));
        for i = floor(3600*1.5/C_Rate):-1:1
            if isnan(V_bat_1(i)) == 1
                V_bat_1(i) = 1;
            end
            if isnan(V_bat_2(i)) == 1
                V_bat_2(i) = 1;
            end
        end

        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
        %% TOLERANZ 
        
        %     tolerance(n) = mean(abs((V_bat_2)./V_bat_1-1)*100);  % Durchschnittliche Abweichung in % für eine Batterie
        tolerance(n) = (flaeche1-flaeche2)/flaeche2*100;
        
    end
    
    tolerance_crate = [tolerance_crate tolerance];
    

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
    %% FEHLERTOLERANZ
    % Plotten der Abweichungen für jede Batterienummer
    x = 1:length(DATA);
    plot(x,tolerance)
    grid on
    xlabel('Batterienummer (Index in Matrix)');
    ylabel('Durchschnittliche Abweichung in %');
    average_tolerance = mean(tolerance);        % durchschnittliche Abweichung in % aller Batterien
end

