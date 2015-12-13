function myclick_2( hObject, eventdata, handles )
%
%   
global GT_move global_x global_y t;

if isempty(GT_move{t}{global_x,3}) % 第一次点的时候，3号cell肯定是空，此时给global_y赋值
    global_y = 3;
end

% 第一个图的操作只出现在前2列
if ~strcmp( get(gco, 'Type'), 'line') || ~isequal( get(gco, 'Color'), [0 1 0])% 必须选中绿线才进行操作而不是背景
%     msgbox('没选中对象，请重选');
else
    if strcmp(get( gco, 'Selected' ), 'off') % 如果线处于未选中状态，则选中它
        set( gco, 'Selected', 'on' );
        % 写入当前位置
        GT_move{t}{global_x,global_y} = [ GT_move{t}{global_x,global_y}, str2double(get(gco, 'DisplayName')) ];
        
    else % 如果已被选中，则取消选中状态
        set( gco, 'Selected', 'off' );
        GT_move{t}{global_x,global_y} = GT_move{t}{global_x,global_y}(1:end-1); % 删去csp中的最后一个
        if isempty( GT_move{t}{global_x,global_y} )
            msgbox('CSP已清空');
        end
        
    end
end
   