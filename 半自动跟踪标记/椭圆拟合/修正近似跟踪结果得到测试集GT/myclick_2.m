function myclick_2( hObject, eventdata, handles)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
global GT_move global_x t;

% �ڶ���ͼ�Ĳ���ֻ�����ڵڶ���,��������
if ~strcmp( get( gco, 'Type' ), 'line') % ����ѡ���߲Ž��в��������Ǳ���
%     msgbox('ûѡ�ж�������ѡ');
else
    if strcmp(get( gco, 'Selected' ), 'off') % ����ߴ���δѡ��״̬����ѡ����
        set( gco, 'Selected', 'on' );
        % ���3��λΪ0����д��3�ţ�����д��4��
        if ~GT_move{t}(global_x,3)
            GT_move{t}(global_x,3) = str2double( get(gco, 'DisplayName') );
        else
            GT_move{t}(global_x,4) = str2double( get(gco, 'DisplayName') );
        end

    else % ����ѱ�ѡ�У���ȡ��ѡ��״̬
        set( gco, 'Selected', 'off' );
%         global_x=global_x-1;
        % ���4��λ��0�������4�ţ��������3��
        if GT_move{t}(global_x,4)
            GT_move{t}(global_x,4) = 0; % �˺�һ�����û�н�������ѡȡ����λ���ϵ������ڴ��ڣ����Ҫ������ 0
        else
            GT_move{t}(global_x,3) = 0;
        end
    end
end
   





end

