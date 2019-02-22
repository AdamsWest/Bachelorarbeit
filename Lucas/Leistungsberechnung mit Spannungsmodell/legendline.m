h1 = plot(1:10);
[hh,icons,plots,txt] = legend({'Line 1'});
p1 = icons(1).Position;
icons(1).Position = [0.3 p1(2) 0];
icons(2).XData = [0.05 0.2];