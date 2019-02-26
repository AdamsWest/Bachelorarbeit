clear
close all 
clc

load('axi_motor_db')

motor_name = axi_motor_db{21,1}; % Motorname eingeben
[K_V, I_0, R_i, m_Mot, S_max, I_max] = ...
Motordata('axi_motor_db',motor_name);
K_V = K_V*2*pi/60;          % Umrechnung in 1/(V*s)
V_max = 4 * S_max;          % maximum voltage (assuming 4 volts per cell)
tau_max = (I_max-I_0)/K_V;      % max torque according to the model of Drela
omega_max = (V_max-R_i*0)*K_V;  % max angular velocity according to the model of Drela

% Grenzen festlegen
Omega = 0:10:round(omega_max/100)*100;
% Omega = Omega*2*pi/60;
tau = 0:0.01:round(tau_max,1);
[X,Y] = meshgrid(Omega,tau);

% Motorkenngrößen
I_Mot = (Y*K_V+I_0);
I_Mot(I_Mot>I_max) = NaN;
U_Mot = (X/K_V+R_i*I_Mot);
U_Mot(U_Mot>V_max) = NaN;
% Eigentlich muesste es auch eine Leistungsgrenze P_max geben.

% Wirkungsgrad
eta_mot = (X.*Y)./(I_Mot.*U_Mot);

% Plot
surf(X*60/(2*pi),Y,eta_mot);
xlabel('\Omega [RPM]');
ylabel('Drehmoment [Nm]');
zlabel('\eta_{Mot}');

fig = gcf;
saveas(gcf,'Motormodell','pdf')