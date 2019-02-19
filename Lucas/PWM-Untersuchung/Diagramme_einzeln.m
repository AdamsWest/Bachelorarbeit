% Restladung über der Höhe
figure(figure_C_Rest_V)
plot(H,C_Rest_V*100,'LineWidth',2);
grid on
hold on
xlabel('Höhe [m]')
ylabel('Restladung der Batterie [%]')
text(6000,90,['Motor = ' motor_name]);
text(6000,80,['Propeller = ' prop_name]);
saveas(gcf,'Flächenflugzeug_C_Rest_V', 'pdf');  


% Drehzahl über der Höhe
figure(figure_omega)
plot(H,Omega/(2*pi)*60,'LineWidth',2)
grid on
hold on
xlabel('Höhe [m]')
ylabel('RPM')
saveas(gcf,'Flächenflugzeug_omega', 'pdf');  


% Motorstrom über der Höhe
figure(figure_I_mot)
plot(H,I_mot,'LineWidth',2)
grid on
hold on
xlabel('Höhe [m]')
ylabel('I_{mot} [A]')
saveas(gcf,'Flächenflugzeug_I_mot', 'pdf');  

% Motorspannung über der Höhe
figure(figure_U_mot)
plot(H,U_mot,'LineWidth',2)
grid on
hold on
xlabel('Höhe [m]')
ylabel('U_{mot} [V]')
saveas(gcf,'Flächenflugzeug_U_mot', 'pdf');  

% Batteriestrom über der Höhe
figure(figure_I_Bat)
plot(H,I_Bat,'LineWidth',2)
grid on
hold on
xlabel('Höhe [m]')
ylabel('I_{Bat} [A]')
saveas(gcf,'Flächenflugzeug_I_Bat', 'pdf');  

% Batteriespannung über der Höhe
figure(figure_U_Bat)
H2 = [0;H];
plot(H2,U_Bat,'LineWidth',2)
grid on
hold on
xlabel('Höhe [m]')
ylabel('U_{Bat} [V]')
saveas(gcf,'Flächenflugzeug_U_Bat', 'pdf'); 

% PWM über der Höhe
figure(figure_PWM)
plot(H,PWM*100,'LineWidth',2)
grid on
hold on
xlabel('Höhe [m]')
ylabel('PWM [%]')
saveas(gcf,'Flächenflugzeug_PWM', 'pdf');  

% Wirkungsgrad
figure(figure_eta)
plot(H,eta_ges*100,'LineWidth',2)
grid on
hold on
xlabel('Höhe [m]')
ylabel('eta_{ges} [%]')
saveas(gcf,'Flächenflugzeug_eta', 'pdf');  

% Steigwinkel Flaechenflugzeug
figure(figure_gamma)
plot(H,gamma_Flaechenflzg,'LineWidth',2)
grid on
hold on
xlabel('Höhe [m]')
ylabel('Bahnneigungswinkel [°]')
saveas(gcf,'Flächenflugzeug_gamma', 'pdf');  

% Geschwindigkeit Flaechenflugzeug
figure(figure_V)
plot(H,V_Flaechenflugzeug,'LineWidth',2)
grid on
hold on
xlabel('Höhe [m]')
ylabel('Fluggeschwindigkeit [m/s]')
saveas(gcf,'Flächenflugzeug_V', 'pdf');