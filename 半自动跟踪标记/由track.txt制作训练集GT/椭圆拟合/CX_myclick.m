function CX_myclick( hObject, eventdata, handles )

global tmp_label2e global_x global_y last_click;

if ~strcmp( get( gco, 'Type' ), 'line') % 必须选中线才进行操作而不是背景
%     msgbox('没选中对象，请重选');
else
    % 2015.8.27 新改动：通过识别当前选中对象（点或线），来进行标记，因此可以正确标记merge事件
    % 2015.10.9 新改动：通过全局变量 last_click 来防止先选线后选点
    % dot_or_line为1时说明刚点击了一个点
    isdot = isequal( get( gco, 'marker' ), '*'); 
    if (global_y==1 && isdot) || (global_y==2 && ~isdot) % 光标在1处且选中点或者光标在2处且选中线，则记录
        set( gco, 'Selected', 'on' );
        last_click = get(gco, 'DisplayName'); % 刚点击完的目标
        tmp_label2e(global_x,global_y) = str2double( last_click );
        % 选中了则光标走起，继续选择下一个
        global_y=global_y+1;
        if global_y==3 % 说明此条记录已完成
            msgbox('完成一次连接！'); 
            global_y=1;
            global_x=global_x+1;
        end
    else % 其他情况：光标1选中线/光标2选中点，都属于删除操作
        if ~strcmp(get(gco, 'DisplayName'), last_click) 
            % 不是点击之前被选中的，则不操作
        else
            set( gco, 'Selected', 'off' );
                if global_y==2 % 光标退后一格，进行重新选取
                    global_y=1;
                    msgbox('删除一条连接！');
                else
                    global_y=2;
                    global_x=global_x-1;   
                    last_click = num2str(tmp_label2e(global_x,1)); % 被选中的也要撤退一格
                end
                if global_x>0 % 第一次就选线时global_x会减到0，需要进行判断
                    tmp_label2e(global_x,global_y) = 0; % 退后一格后若没有进行重新选取，则位置上的数仍在存在，因此要将其置 0
                else
                    error('必须先选择点再选择椭圆！');
                end
        end
    end
    
% 以下内容为旧版，无法标记merge
%     if strcmp(get( gco, 'Selected' ), 'off') % 如果线处于未选中状态，则选中它
%         set( gco, 'Selected', 'on' );
%         tmp_label2e(global_x,global_y) = str2double( get(gco, 'DisplayName') );
%         % 选中了则光标走起，继续选择下一个
%         global_y=global_y+1;
%         if global_y==3
%             global_y=1;
%             global_x=global_x+1;
%         end
%     else % 如果已被选中，则取消选中状态         
%         set( gco, 'Selected', 'off' );
%         if global_y==2 % 光标退后一格，进行重新选取
%             global_y=1;
%         else
%             global_y=2;
%             global_x=global_x-1;
%         end
%         tmp_label2e(global_x,global_y) = 0; % 退后一格后若没有进行重新选取，则位置上的数仍在存在，因此要将其置 0
%     end

end
    
 
end