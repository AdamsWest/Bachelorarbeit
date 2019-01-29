% Motorregression

clc
clear
clf
close

load('axi_motor_db.mat');
format long
[K_V, I_0, R_i, m_Mot, S_max, I_max, ges] = extraction_axi(axi_motor_db);

% Initialisierung

m_Mot2 = m_Mot;
m_Mot2(50) = [];
m_Mot3 = m_Mot2;
m_Mot3(46) = [];
I_max2 = I_max;
I_max2(50) = [];
I_max3 = I_max2;
I_max3(46) = [];
S_max2 = S_max;
S_max2(50) = [];
S_max3 = S_max2;
S_max3(46) = [];
K_V2 = K_V;
K_V2(50) = [];
K_V3 = K_V2;
K_V3(46) = [];

% Faktor b Berechnung

x1 = m_Mot;
y1 = I_max .* (S_max * 3.7);                % P_max
x2 = m_Mot2;
y2 = I_max2 .* (S_max2 * 3.7);              % P_max
b1 = x1\y1;
b2 = x2\y2;
x3 = m_Mot3;
y3 = I_max3 .* (S_max3 * 3.7);
b3 = x3\y3;

yCalc1 = b1*x1;
yCalc2 = b2*x2;
yCalc3 = b3*x3;

figure(1);
hold on
scatter(
plot(x1,yCalc1,'b', x2,yCalc2,'r', x3,yCalc3, 'k')



