function myclick_2( hObject, eventdata, handles )
%
%   
global GT_move global_x global_y t;

if isempty(GT_move{t}{global_x,3}) % ��һ�ε��ʱ��3��cell�϶��ǿգ���ʱ��global_y��ֵ
    global_y = 3;
end

% ��һ��ͼ�Ĳ���ֻ������ǰ2��
if ~strcmp( get(gco, 'Type'), 'line') || ~isequal( get(gco, 'Color'), [0 1 0])% ����ѡ�����߲Ž��в��������Ǳ���
%     msgbox('ûѡ�ж�������ѡ');
else
    if strcmp(get( gco, 'Selected' ), 'off') % ����ߴ���δѡ��״̬����ѡ����
        set( gco, 'Selected', 'on' );
        % д�뵱ǰλ��
        GT_move{t}{global_x,global_y} = [ GT_move{t}{global_x,global_y}, str2double(get(gco, 'DisplayName')) ];
        
    else % ����ѱ�ѡ�У���ȡ��ѡ��״̬
        set( gco, 'Selected', 'off' );
        GT_move{t}{global_x,global_y} = GT_move{t}{global_x,global_y}(1:end-1); % ɾȥcsp�е����һ��
        if isempty( GT_move{t}{global_x,global_y} )
            msgbox('CSP�����');
        end
        
    end
end
   