function press_c_space( hObject, eventdata, handles)
%
%  
global GT_move GT_delete global_x global_y t;

keyPress = get(gcf, 'currentcharacter');

if keyPress=='c'
    % ����'c'���bsp��¼��ϣ��������
    msgbox('csp choose done!');
    global_y = global_y + 1;
end

% ��һ֡�������������ʧ�¼���������Բ���¿ո�����ͻ��¼
if keyPress==' '
    % ���¿ո�����м�¼��ϣ���Ҫ�ж��¼�
    abcd = num2str( ~isemptycell(GT_move{t}(global_x,:)) );
    switch abcd
        case '1  0  0  0'
            msgbox('disappear');
        case '0  0  1  0'
            msgbox('appear');
        case '1  0  1  0'
            msgbox('move');
        case '1  0  1  1'
            msgbox('divide/split');
        case '1  1  1  0'
            msgbox('merge');
        otherwise
            % ��������򽫴�����գ�������ѡȡ
            GT_move{t}(global_x,:) = cell(1,4);
            global_x = global_x - 1;
            msgbox('error! please re choose'); 
    end
    global_x = global_x+1; % ������һ��
    global_y = 1; % ���ظ�
end

if gcf==1 && strcmp(keyPress, 'd')
    % �����Ҫɾ������Բ����d��ֻ���ڵ�һ֡�д�����¼�
    % ���Ϊɾ������Բ���ں������Ƴ���·��ȥ·
    GT_delete{t} = [ GT_delete{t}; GT_move{t}(global_x,1)];
    GT_move{t}{global_x,1} = []; % ��Ҫ��֮ǰ��¼��Ǩ����Ϣ����
    % ��ʾ������delete����
    msgbox('delete');
end



end

