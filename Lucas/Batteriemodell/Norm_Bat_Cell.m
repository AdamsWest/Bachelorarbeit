% Normierung von Elektromodellflug
% clear
close all
clc


load('Elektromodellflug.mat');
Elektromodellflug_norm = Elektromodellflug;


% Initialisierungen
sum = 0;
sum_1 = 0;
sum_2 = 0;
sum_3 = 0;
sum_4 = 0;
sum_5 = 0;
sum_6 = 0;
sum_7 = 0;
sum_8 = 0;


capacity = zeros(length(Elektromodellflug_norm),1);
resistance = zeros(length(Elektromodellflug_norm),1);

for i = 1:length(Elektromodellflug_norm)
    
    % Normierung von Q, Q_nom und Q_exp mit der Kapazität in As
    
    Elektromodellflug_norm{i,3}(1) = Elektromodellflug_norm{i,3}(1) * 1000 / (Elektromodellflug_norm{i,5}*3600);
    Elektromodellflug_norm{i,3}(2) = Elektromodellflug_norm{i,3}(2) * 1000 / (Elektromodellflug_norm{i,5}*3600);
    Elektromodellflug_norm{i,3}(3) = Elektromodellflug_norm{i,3}(3) * 1000 / (Elektromodellflug_norm{i,5}*3600);

    % Normzelle: (arithmetischer Mittelwert)
    
    sum_1 = sum_1 + Elektromodellflug_norm{i,3}(1);
    sum_2 = sum_2 + Elektromodellflug_norm{i,3}(2);
    sum_3 = sum_3 + Elektromodellflug_norm{i,3}(3);
    sum_4 = sum_4 + Elektromodellflug_norm{i,3}(4);
    sum_5 = sum_5 + Elektromodellflug_norm{i,3}(5);
    sum_6 = sum_6 + Elektromodellflug_norm{i,3}(6);
    sum_7 = sum_7 + Elektromodellflug_norm{i,3}(7);
    sum_8 = sum_8 + Elektromodellflug_norm{i,3}(8);
    
    capacity(i) = Elektromodellflug_norm{i,5}/1000;
    resistance(i) = Elektromodellflug_norm{i,3}(8);
    
    % sum = sum + Elektromodellflug_norm{i,5}/(Elektromodellflug_norm{i,4}*3.6);

end


% arithmetischer Mittelwert über alle Batterien
% durchschnittliche Kapazität pro Zelle:
% Cnom = sum / length(Elektromodellflug_norm);
Q = sum_1 / length(Elektromodellflug_norm);
Qnom = sum_2 / length(Elektromodellflug_norm);
Qexp = sum_3 / length(Elektromodellflug_norm);
Vfull = sum_4 / length(Elektromodellflug_norm);
Vexp = sum_5 / length(Elektromodellflug_norm);
Vnom = sum_6 / length(Elektromodellflug_norm);
i = sum_7 / length(Elektromodellflug_norm);
R = sum_8 / length(Elektromodellflug_norm);

% Cnom wurde eingefügt --> beachten
M_A = [1, 1, 0 ; 1, exp(-3), -Q/(Q-Qexp)*(Qexp+i) ; 1, exp(-3*Qnom/Qexp), -Q/(Q-Qnom)*(Qnom + i)];
b = [Vfull + R*i ; Vexp + R*i ; Vnom + R*i];
x = M_A\b;

Eo = x(1);
A = x(2);
K = x(3);


% xdata = log10(capacity);
% ydata = log10(resistance);
% f = fit(xdata,ydata,'poly1');
% plot(f,xdata,ydata)
% hold on
% p1 = -1.129;
% p2 = -1.73;
% 
% figure
% vec = 0.1:0.1:max(capacity);
% func = 10.^(p2+p1*vec);
% % func = 10^(p2)./(10.^(p1 * vec));
% plot(vec,func, capacity, resistance,'rx')

plot(capacity, resistance,'rx');
vec = 0.1:0.01:max(resistance);
func = 1./vec;
plot(vec,func);

