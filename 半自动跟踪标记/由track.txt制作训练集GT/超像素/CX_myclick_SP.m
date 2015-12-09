function CX_myclick_SP( hObject, eventdata, handles )

global tmp_label2e global_x global_y last_click;

if ~strcmp( get( gco, 'Type' ), 'line') % 必须选中线才进行操作而不是背景
%     msgbox('没选中对象，请重选');
else
    % 2015.8.27 新改动：通过识别当前选中对象（点或线），来进行标记，因此可以正确标记merge事件
    % 2015.10.9 新改动：通过全局变量 last_click 来防止先选线后选点
    % dot_or_line为1时说明刚点击了一个点
    isdot = isequal( get( gco, 'marker' ), '*');
    
    % 光标不在1处却选中点，或光标在1处却选中线，都属于误操作
    if global_y~=1 && isdot && ~strcmp(get(gco, 'DisplayName'), last_click) 
        msgbox('确认完成请按空格后再选择下一个点!');
        return
    end
    if global_y==1 && ~isdot
        msgbox('请先选择点!');
        return
    end
        
    if strcmp(get(gco, 'DisplayName'), last_click) && strcmp(get(gco, 'Selected'),'on')
        % 选中了上一次所选的，则进行删除
        set( gco, 'Selected', 'off' );
        global_y = global_y - 1; % 光标退后一格，进行重新选取
        if global_y==1
            msgbox('删除一条连接！');
        end    
    else
        % 不是选中上一次选的，则进行记录
        set( gco, 'Selected', 'on' );
        tmp_label2e(global_x, global_y) = str2double( get(gco, 'DisplayName') ); % 记录
        % 选中了则光标走起，继续选择下一个
        global_y = global_y + 1;
        
        
        if global_y>size(tmp_label2e,2)
            msgbox('选择的Superpixel数目已到上限，请不要再选！');
        end
    end
    
    last_click = get(gco, 'DisplayName'); % 更新last_click

end
    
