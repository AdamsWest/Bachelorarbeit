close all


figure_C_Rest_V = figure;
figure_omega = figure;
figure_I_mot = figure;
figure_U_mot = figure;
figure_I_Bat = figure;
figure_U_Bat = figure;
figure_PWM = figure;
figure_eta = figure;
% figure_gamma = figure;
figure_V = figure;


% Restladung über der Höhe
figure(figure_C_Rest_V)
plot(H,C_Rest_V*100,'LineWidth',2);
% hold on
% plot(H,C_Rest_V_2*100,'LineWidth',2);
grid on
hold on
xlabel('Höhe [m]')
ylabel('Restladung der Batterie [%]')
% text(6000,90,['Motor = ' motor_name]);
% text(6000,80,['Propeller = ' prop_name]);
saveas(gcf,'Quadrocopter_C_Rest_V', 'png');  


% Drehzahl über der Höhe
figure(figure_omega)
plot(H,Omega/(2*pi)*60,'LineWidth',2)
% hold on 
% plot(H,Omega_2/(2*pi)*60,'LineWidth',2)
grid on
hold on
xlabel('Höhe [m]')
ylabel('RPM')
% legend('V_{Kg} = const.', 'V_{Kg} = var.', 'Location','southeast')
saveas(gcf,'Quadrocopter_omega', 'png');  


% Motorstrom über der Höhe
figure(figure_I_mot)
plot(H,I_mot,'LineWidth',2)
grid on
hold on
xlabel('Höhe [m]')
ylabel('I_{mot} [A]')
saveas(gcf,'Quadrocopter_I_mot', 'png');  

% Motorspannung über der Höhe
figure(figure_U_mot)
plot(H,U_mot,'LineWidth',2)
grid on
hold on
xlabel('Höhe [m]')
ylabel('U_{mot} [V]')
saveas(gcf,'Quadrocopter_U_mot', 'png');  

% Batteriestrom über der Höhe
figure(figure_I_Bat)
plot(H,I_Bat,'LineWidth',2)
grid on
hold on
xlabel('Höhe [m]')
ylabel('I_{Bat} [A]')
saveas(gcf,'Quadrocopter_I_Bat', 'png');  

% Batteriespannung über der Höhe
figure(figure_U_Bat)
H2 = [0;H];
plot(H2,U_Bat,'LineWidth',2)
grid on
hold on
xlabel('Höhe [m]')
ylabel('U_{Bat} [V]')
saveas(gcf,'Quadrocopter_U_Bat', 'png'); 

% PWM über der Höhe
figure(figure_PWM)
plot(H,PWM*100,'LineWidth',2)
% hold on 
% plot(H,PWM_2*100,'LineWidth',2)
grid on
hold on
xlabel('Höhe [m]')
ylabel('PWM [%]')
saveas(gcf,'Quadrocopter_PWM', 'png');  

% Wirkungsgrad
figure(figure_eta)
plot(H,eta_ges*100,'LineWidth',2)
grid on
hold on
xlabel('Höhe [m]')
ylabel('eta_{ges} [%]')
saveas(gcf,'Quadrocopter_eta', 'png');  

% % Steigwinkel Flaechenflugzeug
% figure(figure_gamma)
% plot(H,gamma_Flaechenflzg,'LineWidth',2)
% grid on
% hold on
% xlabel('Höhe [m]')
% ylabel('Bahnneigungswinkel [°]')
% saveas(gcf,'Flächenflugzeug_gamma', 'pdf');  
% 
% Geschwindigkeit Flaechenflugzeug
figure(figure_V)
plot(H,V_Kg,'LineWidth',2)
% hold on
% plot(H,V_Kg_2,'LineWidth',2)
grid on
hold on
xlabel('Höhe [m]')
ylabel('Fluggeschwindigkeit [m/s]')
% legend('V_{Kg} = const.', 'V_{Kg} = var.')%,'Location','southeast')
saveas(gcf,'Flächenflugzeug_V', 'png');