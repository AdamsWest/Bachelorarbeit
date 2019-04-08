figure
width = 20;
heigth = 20;

plot(peaks);
ax = gca;
% fig.PaperPositionMode = 'auto';
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_heigth = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_heigth];
set(gcf,'Units','centimeters', 'PaperSize', [width heigth]);
print(gcf, 'Fig.pdf', '-dpdf', '-fillpage')