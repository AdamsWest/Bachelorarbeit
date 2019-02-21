figure_geschw = figure;
figure(figure_geschw);
subplot(321), stairs(H,vi,'LineWidth',1), grid, title('Geschwindigkeiten'), xlabel('Höhe [m]'),ylabel('[m/s]')
    hold on
    stairs(H,mu_z,'LineWidth',1)
    legend('induz. Geschw. v_i','\mu_z','Location','southeast');
    hold off
subplot(322), stairs(H,Thrust,'LineWidth',1), grid, title('Schub'), xlabel('Höhe [m]'),ylabel('Schub [N]')
subplot(323), stairs(H,Omega/(2*pi)*60,'LineWidth',1), grid, title('Drehzahl'), xlabel('Höhe [m]'),ylabel('Drehzahl [RPM]')
    hold on 
    stairs(H,Omega/(2*pi)*60./Uebersetzung,'LineWidth',1)
    legend('Motordrehzahl', 'Propellerdrehzahl','Location','southeast')
subplot(324), stairs(H,tau,'LineWidth',1), grid, title('Drehmoment'), xlabel('Höhe [m]'),ylabel('Drehmoment [Nm]')
    hold on 
    stairs(H,tau .* Uebersetzung,'LineWidth',1), legend('Motordrehmoment','Propellerdrehmoment','Location','northeast');
subplot(325), stairs(H,Omega.*tau,'LineWidth',1), grid, title('Gegenüberstellung der Leistungen'), xlabel('Höhe [m]'), ylabel('Leistung [W]');
    hold on 
    stairs(H,Thrust.*(mu_z+vi),'LineWidth',1), legend('Wellenleistung (M*\omega)','Strahlleistung (T*(v_i+\mu_Z)','Location','southwest');

    ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'MySavedFile','-dpdf')