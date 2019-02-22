x = 0:0.5:10;
figure; hold on;
ph(1) = plot(x,sin(x), 'Marker', 'o');
ph(2) = plot(x,cos(x), 'Marker', 's');

ax = gca;
ax.Box = 'on';

bx = axes('position', [ax.Position(1:2), 0.3, 0.4], 'Box', 'on', 'XTick', [], 'YTick', []);%, 'Color', [0.8549,0.7020,1.0000]);
cx = axes('position', [ax.Position(1:2), 0.3, 0.4], 'Box', 'on', 'XTick', [], 'YTick', []);%, 'Color', [0.8549,0.5020,1.0000]);

[legb, objsb] = legend(bx, ph, {'sin', 'cos'}, 'Location', 'SouthWest');
[legc, objsc] = legend(cx, ph, {'sin', 'cos'}, 'Location', 'SouthWest');

line_start_end = [0.01, 0.3];
line_text_step = 0.05;

legendshrink(ax.Position(1:2), legb, objsb, bx, line_start_end, line_text_step);
legendshrink(ax.Position(1:2), legc, objsc, cx, line_start_end, line_text_step);

% you need only to adjust cx' position.
cx.Position(1:2) = [ax.Position(1)+ax.Position(3)-cx.Position(3), ax.Position(2)];
legc.Position(1:2) = cx.Position(1:2);
