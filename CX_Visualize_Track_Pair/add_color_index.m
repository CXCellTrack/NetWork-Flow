function e = add_color_index( e )

% ����һ���쳣�����
if ~isfield(e,'color_index')
    e.color_index = randi(64);
end