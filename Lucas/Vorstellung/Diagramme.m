% Darstellung der Ergenisse in Diagrammen
figure_C_Rest_V = figure;
figure_omega = figure;
figure_M_tip = figure;
figure_eta_ges = figure;



% Restladung über der Höhe
figure(figure_C_Rest_V)
plot(H,C_Rest_V_1*100,'LineWidth',2);
grid on
hold on
plot(H,C_Rest_V_2*100,'r-','LineWidth',2);
xlabel('Höhe [m]')
ylabel('Restladung der Batterie [%]')


% Drehzahl über der Höhe
figure(figure_omega)
plot(H,Omega_1/(2*pi)*60,'LineWidth',2)
grid on
hold on
plot(H,Omega_2/(2*pi)*60,'r-','LineWidth',2);
xlabel('Höhe [m]')
ylabel('RPM')


% Drehzahl über der Höhe
figure(figure_M_tip)
plot(H,M_tip_1,'LineWidth',2)
grid on
hold on
plot(H,M_tip_2,'r-','LineWidth',2);
xlabel('Höhe [m]')
ylabel('M_{tip}')


% Wirkungsgrad
figure(figure_eta_ges)
plot(H,eta_ges_1*100,'LineWidth',2)
grid on
hold on
plot(H,eta_ges_2*100,'r-','LineWidth',2);
xlabel('Höhe [m]')
ylabel('eta_{ges} [%]')

%% Datei abspeichern
%ImageSizeX = 40;
%ImageSizeY = 30;
figure(figure_C_Rest_V)
%set(gcf,'PaperUnits','centimeters', 'PaperPosition', [0 0 ImageSizeX ImageSizeY]); 
%set(gcf,'Units','centimeters', 'PaperSize', [ImageSizeX ImageSizeY]); 
saveas(gcf,'C_Rest_V', 'jpg'); 
figure(figure_omega)
saveas(gcf,'omega', 'jpg');
%figure(figure_I_mot)
%saveas(gcf,'I_mot', 'jpg');
%figure(figure_U_mot)
%saveas(gcf,'U_mot', 'jpg');
%figure(figure_I_Bat)
%saveas(gcf,'I_Bat', 'jpg');
%figure(figure_PWM)
%saveas(gcf,'PWM', 'jpg');
figure(figure_M_tip)
saveas(gcf,'M_tip', 'jpg');
figure(figure_eta_ges)
saveas(gcf,'eta_ges', 'jpg');