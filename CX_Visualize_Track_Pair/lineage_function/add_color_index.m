function e = add_color_index( e )

% 这是一种异常情况！
if ~isfield(e,'color_index')
    e.color_index = randi(64);
end