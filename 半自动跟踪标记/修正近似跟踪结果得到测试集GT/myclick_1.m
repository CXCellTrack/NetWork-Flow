function myclick_1( hObject, eventdata, handles)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
global GT_move global_x t;

% ��һ��ͼ�Ĳ���ֻ�����ڵ�һ��
if ~strcmp( get(gco, 'Type'), 'line') || ~isequal( get(gco, 'Color'), [0 1 0])% ����ѡ�����߲Ž��в��������Ǳ���
%     msgbox('ûѡ�ж�������ѡ');
else
    if strcmp(get( gco, 'Selected' ), 'off') % ����ߴ���δѡ��״̬����ѡ����
        set( gco, 'Selected', 'on' );
        % ���1��λΪ0����д��1�ţ�����д��2��
        if ~GT_move{t}(global_x,1)
            GT_move{t}(global_x,1) = str2double( get(gco, 'DisplayName') );
        else
            GT_move{t}(global_x,2) = str2double( get(gco, 'DisplayName') );
        end
        
    else % ����ѱ�ѡ�У���ȡ��ѡ��״̬
        set( gco, 'Selected', 'off' );
        % ���2��λ��0�������2�ţ��������1��
        if GT_move{t}(global_x,2)
            GT_move{t}(global_x,2) = 0; % �˺�һ�����û�н�������ѡȡ����λ���ϵ������ڴ��ڣ����Ҫ������ 0
        else
            GT_move{t}(global_x,1) = 0;
        end
    end
end
   





end

