function pressBlank( hObject, eventdata, handles)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
global GT_move GT_delete global_x t;

keyPress = get(gcf, 'currentcharacter');


% ��һ֡�������������ʧ�¼���������Բ���¿ո�����ͻ��¼
if keyPress==' '

    % ���¿ո�����м�¼��ϣ���Ҫ�ж��¼�
    abcd = num2str(GT_move{t}(global_x,:)~=0);
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
            GT_move{t}(global_x,:) = [0 0 0 0];
            global_x = global_x - 1;
            msgbox('error! please re choose'); 
    end
    global_x=global_x+1;
end

if gcf==1 && strcmp(keyPress, 'd')
    % �����Ҫɾ������Բ����d��ֻ���ڵ�һ֡�д�����¼�
    % ���Ϊɾ������Բ���ں������Ƴ���·��ȥ·
    GT_delete{t} = [ GT_delete{t}; GT_move{t}(global_x,1)];
    GT_move{t}(global_x,1) = 0; % ��Ҫ��֮ǰ��¼��Ǩ����Ϣ����
    % ��ʾ������delete����
    msgbox('delete');
end



end

