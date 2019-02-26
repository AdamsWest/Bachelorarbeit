% Propellerkennfeldplot
clear
close all
clc
load('DATA_APC');

[RPM, V, T, P, Tau] = Propeller_map(DATA_APC,'10x3');

surf(V,RPM,T)
xlabel('V [m/s]')
ylabel('\Omega [RPM]')
zlabel('Schub [N]')
zlim([min(min(T)) 100])
saveas(gcf,'Propellerkennfeld','pdf')
% [RPM,V,T,P]=Propeller_map(DATA_APC,'10x3');     % Propellerkennfeld erzeugen 
% eta = T.*repmat(V,length(T(:,1)),1)./P;         % Wirkungsgrad = T*V / P (ohne V_i) berechnen 
% figure 
% surf(V,RPM,eta) 
% xlabel('V, m/s') 
% ylabel('\omega, RPM') 
% zlabel('\eta, -') 
% zlim([0 1]) 