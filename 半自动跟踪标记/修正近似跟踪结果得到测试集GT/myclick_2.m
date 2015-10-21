function myclick_2( hObject, eventdata, handles)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
global GT_move global_x t;

% 第二个图的操作只出现在第二列,并负责换行
if ~strcmp( get( gco, 'Type' ), 'line') % 必须选中线才进行操作而不是背景
%     msgbox('没选中对象，请重选');
else
    if strcmp(get( gco, 'Selected' ), 'off') % 如果线处于未选中状态，则选中它
        set( gco, 'Selected', 'on' );
        % 如果3号位为0，则写入3号，否则写入4号
        if ~GT_move{t}(global_x,3)
            GT_move{t}(global_x,3) = str2double( get(gco, 'DisplayName') );
        else
            GT_move{t}(global_x,4) = str2double( get(gco, 'DisplayName') );
        end

    else % 如果已被选中，则取消选中状态
        set( gco, 'Selected', 'off' );
%         global_x=global_x-1;
        % 如果4号位非0，则清空4号，否则清空3号
        if GT_move{t}(global_x,4)
            GT_move{t}(global_x,4) = 0; % 退后一格后若没有进行重新选取，则位置上的数仍在存在，因此要将其置 0
        else
            GT_move{t}(global_x,3) = 0;
        end
    end
end
   





end

