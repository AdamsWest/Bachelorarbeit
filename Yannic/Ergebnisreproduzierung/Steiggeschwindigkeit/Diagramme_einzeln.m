% Restladung �ber der H�he
figure(figure_C_Rest_V)
plot(H,C_Rest_V*100,'LineWidth',2);
grid on
hold on
xlabel('H�he [m]')
ylabel('Restladung der Batterie [%]')
text(6000,90,['Motor = ' motor_name]);
text(6000,80,['Propeller = ' prop_name]);
saveas(gcf,'Fl�chenflugzeug_C_Rest_V', 'pdf');  


% Drehzahl �ber der H�he
figure(figure_omega)
plot(H,Omega/(2*pi)*60,'LineWidth',2)
grid on
hold on
xlabel('H�he [m]')
ylabel('RPM')
saveas(gcf,'Fl�chenflugzeug_omega', 'pdf');  


% Motorstrom �ber der H�he
figure(figure_I_mot)
plot(H,I_mot,'LineWidth',2)
grid on
hold on
xlabel('H�he [m]')
ylabel('I_{mot} [A]')
saveas(gcf,'Fl�chenflugzeug_I_mot', 'pdf');  

% Motorspannung �ber der H�he
figure(figure_U_mot)
plot(H,U_mot,'LineWidth',2)
grid on
hold on
xlabel('H�he [m]')
ylabel('U_{mot} [V]')
saveas(gcf,'Fl�chenflugzeug_U_mot', 'pdf');  

% Batteriestrom �ber der H�he
figure(figure_I_Bat)
plot(H,I_Bat,'LineWidth',2)
grid on
hold on
xlabel('H�he [m]')
ylabel('I_{Bat} [A]')
saveas(gcf,'Fl�chenflugzeug_I_Bat', 'pdf');  

% Batteriespannung �ber der H�he
figure(figure_U_Bat)
H2 = [0;H];
plot(H2,U_Bat,'LineWidth',2)
grid on
hold on
xlabel('H�he [m]')
ylabel('U_{Bat} [V]')
saveas(gcf,'Fl�chenflugzeug_U_Bat', 'pdf'); 

% PWM �ber der H�he
figure(figure_PWM)
plot(H,PWM*100,'LineWidth',2)
grid on
hold on
xlabel('H�he [m]')
ylabel('PWM [%]')
saveas(gcf,'Fl�chenflugzeug_PWM', 'pdf');  

% Wirkungsgrad
figure(figure_eta)
plot(H,eta_ges*100,'LineWidth',2)
grid on
hold on
xlabel('H�he [m]')
ylabel('eta_{ges} [%]')
saveas(gcf,'Fl�chenflugzeug_eta', 'pdf');  

% Steigwinkel Flaechenflugzeug
figure(figure_gamma)
plot(H,gamma_Flaechenflzg,'LineWidth',2)
grid on
hold on
xlabel('H�he [m]')
ylabel('Bahnneigungswinkel [�]')
saveas(gcf,'Fl�chenflugzeug_gamma', 'pdf');  

% Geschwindigkeit Flaechenflugzeug
figure(figure_V)
plot(H,V_Flaechenflugzeug,'LineWidth',2)
grid on
hold on
xlabel('H�he [m]')
ylabel('Fluggeschwindigkeit [m/s]')
saveas(gcf,'Fl�chenflugzeug_V', 'pdf');