close all
figure_vpp = figure;

figure(figure_vpp);

subplot(621), stairs(H,C_Rest_V*100,'LineWidth',1), grid, title('Restladung'), hold on
stairs(H,C_Rest_V_vpp*100,'LineWidth',1), xlabel('Höhe [m]'),ylabel('C_{Bat,Rest} [%]')
subplot(622), stairs(H,Omega/(2*pi)*60,'LineWidth',1), grid, title('Drehzahl'), hold on
stairs(H,Omega_vpp/(2*pi)*60,'LineWidth',1), xlabel('Höhe [m]'),ylabel('Drehzahl [RPM]')
subplot(623), stairs(H,I_mot,'LineWidth',1), grid, title('Motorstrom'), hold on
stairs(H,I_mot_vpp,'LineWidth',1), xlabel('Höhe [m]'),ylabel('I_{Mot} [A]')
subplot(624), stairs(H,U_mot,'LineWidth',1), grid,  title('Motorspannung'), hold on
stairs(H,U_mot_vpp,'LineWidth',1), xlabel('Höhe [m]'),ylabel('U_{mot} [V]')
subplot(625), stairs(H,I_Bat,'LineWidth',1), grid, title('Batteriestrom'), hold on
stairs(H,I_Bat_vpp,'LineWidth',1), xlabel('Höhe [m]'),ylabel('I_{Bat} [A]')
H2 = [0;H];
subplot(626), stairs(H2,U_Bat,'LineWidth',1), grid, title('Batteriespannung'), hold on
stairs(H2,U_Bat_vpp,'LineWidth',1), xlabel('Höhe [m]'),ylabel('U_{Bat} [V]')
subplot(627), stairs(H,PWM*100,'LineWidth',1), grid, title('Pulsweitenmodulation'), hold on
stairs(H,PWM_vpp*100,'LineWidth',1), xlabel('Höhe [m]'),ylabel('PWM [%]')
subplot(628), stairs(H,eta_ges*100,'LineWidth',1), grid, title('Gesamtwirkungsgrad'), hold on
stairs(H,eta_ges_vpp*100,'LineWidth',1), xlabel('Höhe [m]'),ylabel('\eta_{ges} [%]')
subplot(629), stairs(H,V_Kg,'LineWidth',1), title('Bahngeschwindigkeit'), grid, hold on
stairs(H,V_Kg_vpp,'LineWidth',1), xlabel('Höhe [m]'),ylabel('V_{Kg} [m/s]')
subplot(6,2,10), stairs(H2,t_Flug,'LineWidth',1), grid, hold on
stairs(H2,t_Flug_vpp,'LineWidth',1), title('Flugzeit'), grid, xlabel('Höhe [m]'),ylabel('t_{Flug} [s]')
pitch = zeros(length(pitch_vpp),1);
for i = 1:length(pitch_vpp)
    if i <= 147
        pitch(i) = 3;
    else
        pitch(i) = NaN;
    end
end
subplot(6,2,11), stairs(H,pitch,'LineWidth',1), hold on
stairs(H,pitch_vpp,'LineWidth',1), title('Pitch'), grid, xlabel('Höhe [m]'),ylabel('Pitch [in]')
lgd = legend('fixed pitch','variable pitch','Location','bestoutside'); title(lgd,'Propeller')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(figure_vpp)
PaperSizeX = 21;
PaperSizeY = 29.7;

fig = gcf;
set(gcf,'PaperUnits','centimeters', 'PaperPosition', [-1.75 -2.1 24.65 34.45]);%[-1.75 -2.6 24.65 34.45]);
set(gcf,'Units','centimeters', 'PaperSize', [PaperSizeX PaperSizeY]);
saveas(gcf,'Verstellpropeller', 'pdf');

