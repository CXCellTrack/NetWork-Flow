function pressBlank_SP( hObject, eventdata, handles)
%
%
global global_x global_y;

keyPress = get(gcf, 'currentcharacter');

% ���¿ո������һ����¼�����л���
if keyPress==' '
    % ���ֻѡ��һ�����������¼
    if global_y==2
        msgbox('��û��ѡ��õ��Ӧ��Superpixel��'); 
        return
    end
    % ���У�����Ƶ���һλ
    global_x = global_x + 1;
    global_y = 1;
    msgbox('���һ�����ӣ�'); 
end


end

