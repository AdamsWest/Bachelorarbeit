% Erstellt lineare Regression durch die Motorpunkte
clc
clear
clf
close all

load('axi_motor_db.mat');
format long
[K_V, I_0, R_i, m_Mot, S_max, I_max, ges] = extraction_axi(axi_motor_db);

m_Mot2 = m_Mot;
m_Mot2(50) = [];
I_max2 = I_max;
I_max2(50) = [];
S_max2 = S_max;
S_max2(50) = [];




x1 = m_Mot;
y1 = I_max .* (S_max * 3.7);
x2 = m_Mot2;
y2 = I_max2 .* (S_max2 * 3.7);
b1 = x1\y1;
b2 = x2\y2;
yCalc1 = b1*x1;
yCalc2 = b2*x2;
scatter(x1,y1)
hold on 
figure(1);
plot(x1,yCalc1,'b', x2,yCalc2,'r')

m_Mot3 = m_Mot2;
I_max3 = I_max2;
S_max3 = S_max2;
m_Mot3(46) = [];
I_max3(46) = [];
S_max3(46) = [];


    
    

plot(x3,yCalc5,'k')


X1 = [ones(length(x1),1),x1];
b3 = X1\y1;
yCalc3 = X1 * b3;
X2 = [ones(length(x2),1),x2];
b4 = X2\y2;
yCalc4 = X2 * b4;

x3 = m_Mot3;
y3 = I_max3 .* (S_max3 * 3.7);
b5 = x3\y3;
yCalc5 = b5*x3;

plot(x1,yCalc3,'--b', x2,yCalc4,'--r')
xlabel('Motorgewicht in g');
ylabel('max. Leistung in W');
legend('Data', 'lin. Regression mit Ausreiﬂer und ohne y_o', 'lin. Regression ohne Ausreiﬂer und ohne y_o', ...
    'lin. Regression mit Ausreiﬂer und mit y_o', 'lin. Regression ohne Ausreiﬂer und mit y_o', ...
    'lin. Regression ohne zwei Ausreiﬂer und ohne y_o');

hold off

K_M1 = 1 / K_V * (30 /pi);
K_V2 = K_V;
K_V2(50) = [];
K_V3 = K_V2;
K_V2(46) = [];



figure(2);
%K_M2 = 1 / K_V2 * (30 /pi);
%K_M3 = 1 / K_V3 * (30 /pi);

plot3(x1, y1, K_M1, 'bx')
xlabel('Motormasse in g');
ylabel('max.Leistung in W');
zlabel('Motorkonstante in Nm/A');


