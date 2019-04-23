close all
figure_ueber = figure;

figure(figure_ueber);

subplot(621), stairs(H,C_Rest_V*100,'LineWidth',1), grid, title('Restladung'), hold on
stairs(H,C_Rest_V_g*100,'LineWidth',1), xlabel('Höhe [m]'),ylabel('C_{Bat,Rest} [%]')
subplot(622), stairs(H,Omega/(2*pi)*60,'LineWidth',1), grid, title('Propellerdrehzahl'), hold on
stairs(H,Omega_g/(2*pi)*60./Uebersetzung,'LineWidth',1), xlabel('Höhe [m]'),ylabel('\Omega [RPM]')
subplot(623), stairs(H,I_mot,'LineWidth',1), grid, title('Motorstrom'), hold on
stairs(H,I_mot_g,'LineWidth',1), xlabel('Höhe [m]'),ylabel('I_{Mot} [A]')
subplot(624), stairs(H,U_mot,'LineWidth',1), grid,  title('Motorspannung'), hold on
stairs(H,U_mot_g,'LineWidth',1), xlabel('Höhe [m]'),ylabel('U_{mot} [V]')
subplot(625), stairs(H,I_Bat,'LineWidth',1), grid, title('Batteriestrom'), hold on
stairs(H,I_Bat_g,'LineWidth',1), xlabel('Höhe [m]'),ylabel('I_{Bat} [A]')
H2 = [0;H];
subplot(626), stairs(H2,U_Bat,'LineWidth',1), grid, title('Batteriespannung'), hold on
stairs(H2,U_Bat_g,'LineWidth',1), xlabel('Höhe [m]'),ylabel('U_{Bat} [V]')
subplot(627), stairs(H,PWM*100,'LineWidth',1), grid, title('Pulsweitenmodulation'), hold on
stairs(H,PWM_g*100,'LineWidth',1), xlabel('Höhe [m]'),ylabel('PWM [%]')
subplot(628), stairs(H,eta_ges*100,'LineWidth',1), grid, title('Gesamtwirkungsgrad'), hold on
stairs(H,eta_ges_g*100,'LineWidth',1), xlabel('Höhe [m]'),ylabel('\eta_{ges} [%]')
subplot(629), stairs(H,V_Kg,'LineWidth',1), title('Bahngeschwindigkeit'), grid, hold on
stairs(H,V_Kg_g,'LineWidth',1), xlabel('Höhe [m]'),ylabel('V_{Kg} [m/s]')
subplot(6,2,10), stairs(H2,t_Flug,'LineWidth',1), grid, hold on
stairs(H2,t_Flug_g,'LineWidth',1), title('Flugzeit'), grid, xlabel('Höhe [m]'),ylabel('t_{Flug} [s]')
ueber = zeros(length(Uebersetzung),1);
for i = 1:length(ueber)
    if i <= 166
        ueber(i) = 1;
    else
        ueber(i) = NaN;
    end
end
subplot(6,2,11), stairs(H,ueber,'LineWidth',1), hold on
stairs(H,Uebersetzung,'LineWidth',1), title('Übersetzung'), grid, xlabel('Höhe [m]'),ylabel('Übersetzung i [-]')
lgd = legend('ohne','mit','Location','bestoutside'); title(lgd,'Getriebe')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(figure_ueber)
% PaperSizeX = 21;
% PaperSizeY = 29.7;
% 
% fig = gcf;
% set(gcf,'PaperUnits','centimeters', 'PaperPosition', [-1.75 -2.1 24.65 34.45]);%[-1.75 -2.6 24.65 34.45]);
% set(gcf,'Units','centimeters', 'PaperSize', [PaperSizeX PaperSizeY]);
% saveas(gcf,'Getriebe', 'pdf');

%%
figure_dud = figure;
figure(figure_dud);
subplot(211), stairs(H,Omega/(2*pi)*60,'LineWidth',1), grid, title('Drehzahl'), xlabel('Höhe [m]'),ylabel('\Omega [RPM]')
    hold on 
    stairs(H,Omega/(2*pi)*60./Uebersetzung,'LineWidth',1)
    legend('Motordrehzahl', 'Propellerdrehzahl','Location','southeast')
subplot(212), stairs(H,tau,'LineWidth',1), grid, title('Drehmoment'), xlabel('Höhe [m]'),ylabel('Drehmoment [Nm]')
    hold on 
    stairs(H,tau .* Uebersetzung,'LineWidth',1), legend('Motordrehmoment','Propellerdrehmoment','Location','northeast');

PaperSizeX = 14.8;
PaperSizeY = 21;

fig = gcf;
set(gcf,'PaperUnits','centimeters', 'PaperPosition', [-1.75 -2.1 24.65 34.45]);%[-1.75 -2.6 24.65 34.45]);
set(gcf,'Units','centimeters', 'PaperSize', [PaperSizeX PaperSizeY]);
saveas(gcf,'DuD_KV', 'pdf');

figure_C_Rest_V = figure;
figure(figure_C_Rest_V)
stairs(H,C_Rest_V*100,'LineWidth',1), grid, title('Restladung'), hold on
stairs(H,C_Rest_V_g*100,'LineWidth',1), xlabel('Höhe [m]'),ylabel('C_{Bat,Rest} [%]')
lgd = legend('ohne','mit','Location','northeast'); title(lgd,'Getriebe')
saveas(gcf,'Getriebe','png')

figure_uebersetzung = figure;
figure(figure_uebersetzung)
stairs(H,ueber,'LineWidth',1), hold on, stairs(H,Uebersetzung,'LineWidth',1), 
title('Übersetzung'), grid, xlabel('Höhe [m]'),ylabel('Übersetzung i [-]')
lgd = legend('ohne','mit','Location','northeast'); title(lgd,'Getriebe')
saveas(gcf,'Getriebe_ueber','png')