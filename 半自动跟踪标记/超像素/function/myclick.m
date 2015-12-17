function myclick( hObject, eventdata, handles )
%
%   
global GT_move global_x global_y t;

if gcf==1 && isempty(GT_move{t}{global_x,1}) % 图1第一次点击的时候，1号cell肯定是空，此时给global_y赋值
    global_y = 1;
end
if gcf==2 && isempty(GT_move{t}{global_x,3}) % 图2第一次点击的时候，3号cell肯定是空，此时给global_y赋值
    global_y = 3;
end


% 第一个图的操作只出现在前2列
if strcmp( get(gco, 'Type'), 'line') && strcmp( get(gco, 'linestyle'), '-') % 必须选中线才进行操作而不是*
    
    % 如果线处于选中状态，且位于当前cell内，则取消选中
    if strcmp(get( gco, 'Selected' ), 'on') && ~isempty(GT_move{t}{global_x,global_y}) && GT_move{t}{global_x,global_y}(end)==str2double(get(gco,'DisplayName'))
        
        set( gco, 'Selected', 'off' );
        GT_move{t}{global_x,global_y} = GT_move{t}{global_x,global_y}(1:end-1); % 删去csp中的最后一个
        if isempty( GT_move{t}{global_x,global_y} )
            msgbox('CSP已清空');
        end
        
    else % 如果未选中，则选中并记录
        set( gco, 'Selected', 'on' );
        % 写入当前位置
        if all(GT_move{t}{global_x,global_y}~=str2double(get(gco, 'DisplayName')))
            GT_move{t}{global_x,global_y} = [ GT_move{t}{global_x,global_y}, str2double(get(gco, 'DisplayName')) ];
        end
        
    end
    
end
   





end

