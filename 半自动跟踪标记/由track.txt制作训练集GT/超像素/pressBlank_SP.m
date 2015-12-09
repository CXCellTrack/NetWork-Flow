function pressBlank_SP( hObject, eventdata, handles)
%
%
global global_x global_y;

keyPress = get(gcf, 'currentcharacter');

% 按下空格则完成一条记录，进行换行
if keyPress==' '
    % 如果只选了一个点则不允许记录
    if global_y==2
        msgbox('还没有选择该点对应的Superpixel！'); 
        return
    end
    % 换行，光标移到第一位
    global_x = global_x + 1;
    global_y = 1;
    msgbox('完成一次连接！'); 
end


end

