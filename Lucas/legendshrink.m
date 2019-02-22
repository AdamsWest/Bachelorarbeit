function legendshrink(leg_pos_xy, leg, objs, bx, line_start_end, line_text_step)

% leg_post_xy = ax.Position(1:2);
% [leg, objs] = legend(bx, line_handles, text_cell);
% line_start_end = [0.01, 0.4];
% line_text_step = 0.01;

% for each line, text object, adjust their position in legend
for i = 1:numel(objs)
  if strcmp(get(objs(i), 'Type'), 'line')
    % line object
    if 2 == numel(get(objs(i), 'XData')) % line 
      set(objs(i), 'XData', line_start_end);
    else % marker on line
      set(objs(i), 'XData', sum(line_start_end)/2);
    end
  else
    %text object
    text_pos = get(objs(i), 'Position');
    text_pos(1) = line_start_end(2) + line_text_step;
    set(objs(i), 'Position', text_pos);
  end
end

% get minimum possible width and height
legend_width_no_right = 0;
for i = 1:numel(objs)
  % legend margin left
  if strcmp(get(objs(i), 'Type'), 'line')
    if numel(get(objs(i), 'XData')) == 2
      leg_margin_x = get(objs(i), 'XData');
      leg_margin_x = leg_margin_x(1)*leg.Position(3);
    end
  else
    cur_right = get(objs(i), 'Extent');
    cur_right = (cur_right(1)+cur_right(3))*leg.Position(3);
    if cur_right > legend_width_no_right
      legend_width_no_right = cur_right;
    end
  end
end
legend_width  = legend_width_no_right + leg_margin_x;
legend_height = leg.Position(4);

bx.Position = [leg_pos_xy, legend_width, legend_height];
leg.Position(1:2) = bx.Position(1:2);
leg.Box = 'off';