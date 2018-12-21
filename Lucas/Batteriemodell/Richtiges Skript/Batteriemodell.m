%% Batteriemodell
% 
% Vergleich zweier Batterien (Originale Zelle zur normierten Zelle)
% 
% 



clear
close all
clc

%% Allgemeine Parameter
load('Elektromodellflug.mat');



id_bat = 30;    % Anmerkung: id_bat 33 zu geringe Spannung
PWM = 0.80;
eta_PWM = 0.7;
I_mot = 8;
n_Prop = 4;
C_Rate = 20;
  
%% ORIGINALZELLE


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
V_bat_1 = zeros(3600*1.5/C_Rate,1);
T = 2;
bar = zeros(3600*1.5/C_Rate,1);

% I_bat = PWM * I_mot / eta_PWM * n_Prop;
    I_bat = C_Rate*C;

control = 0;

for n = 1:3600*1.5/C_Rate
    
    
    i_int = I_bat*step_dt + i_int;   % integral of the current
    
    
    % calculating the battery voltage (of one cell)
    V_bat_1(n) = Eo - R*I_bat - K * Q / (Q - i_int) * (i_int + I_bat*0) + A * exp(-B*i_int);
    V_bat_1(n) = N_el * V_bat_1(n);   % the battery voltage of all cells
    
    
    if V_bat_1 < (3.1 * N_el)
        V_bat_1 = NaN;
        % break
    end
    
    
        
    if V_bat_1(n) < 3.1*N_el
        control = 1;
    end
    if control == 1
        V_bat_1(n) = NaN;
    end
    
    
    
    bar(n) = 3.1*N_el;
end
% Kurve der Batterie plotten
x = 1:3600*1.5/C_Rate;
plot(x,V_bat_1)
hold on
plot(x,bar)
hold on



%% NORMZELLE ERZEUGEN

run('Norm_Bat_Cell')
Cnom = Elektromodellflug{id_bat,5}/1000;
B = 3/Qexp;
R = 0.05;

%%
% I_bat = PWM * I_mot / eta_PWM * n_Prop;                         % Batteriestrom
I_bat = C_Rate*C;
I_bat = I_bat/(Cnom);                                           % Normierung des Batteriestroms
i_int = 0;
V_bat_2 = zeros(3600*1.5/C_Rate,1);

control = 0;


for n = 1:3600*1.5/C_Rate
    
    
    i_int = I_bat*step_dt/3600 + i_int;   % <-- HIER ist teilen durch 3600 notwendig, um auf eine Stunde zu skalieren, integral of the current
    
    
    % calculating the battery voltage (of one cell)
    V_bat_2(n) = Eo - R*I_bat - K * Q / (Q - i_int) * (i_int + I_bat*0) + A * exp(-B*i_int);
    % the battery voltage of all cells
    V_bat_2(n) = N_el * V_bat_2(n);  % <-- KORREKTURFAKTOR VON DURCHSCHNITTLICH 5%
    
    
    if V_bat_2 < 3.1 * N_el
        V_bat_2 = NaN;
        break
    end
    
    
    if V_bat_2(n) < 3.1*N_el
        control = 1;
    end
    if control == 1
        V_bat_2(n) = NaN;
    end
    
    
end


% Berechnung der Toleranz nur bis zum ersten Mal, wenn V_bat unter
% V_bat < 3.1 * N_el, ansonsten alle Ergebnisse entfernen
flaeche1 = trapz(V_bat_1(~isnan(V_bat_1)));
flaeche2 = trapz(V_bat_2(~isnan(V_bat_2)));
% for n = floor(3600*1.5/C_Rate):-1:1
%     if isnan(V_bat_1(n)) == 1
%         V_bat_1(n) = 1;
%     end
%     if isnan(V_bat_2(n)) == 1
%         V_bat_2(n) = 1;
%     end
% end

tolerance = (flaeche1-flaeche2)/flaeche2*100;




x = 1:3600*1.5/C_Rate;
plot(x,V_bat_2,'r-')
xlim([0 3600*1.1/C_Rate])
ylim([0 N_el*4.5])
legend('Originaldaten','Grenze','Normierungsmodell')




%% Fehlertoleranz
% Plotten der Toleranz für 
% figure
% a = abs((V_bat_2)./V_bat_1-1)*100;
% plot(x,a)
