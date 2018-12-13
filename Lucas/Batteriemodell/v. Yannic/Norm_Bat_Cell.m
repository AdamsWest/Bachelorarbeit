% Normierung von Elektromodellflug
% clear
% close all
clc


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
% plot(capacity,resistance,'rx')
% hold on
% 
% point_x = 1;
% point_y = 0.015;
% plot(point_x, point_y,'bx')

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




