function myclick( hObject, eventdata, handles )
%
%   
global GT_move global_x global_y t;

if gcf==1 && isempty(GT_move{t}{global_x,1}) % ͼ1��һ�ε����ʱ��1��cell�϶��ǿգ���ʱ��global_y��ֵ
    global_y = 1;
end
if gcf==2 && isempty(GT_move{t}{global_x,3}) % ͼ2��һ�ε����ʱ��3��cell�϶��ǿգ���ʱ��global_y��ֵ
    global_y = 3;
end


% ��һ��ͼ�Ĳ���ֻ������ǰ2��
if strcmp( get(gco, 'Type'), 'line') && strcmp( get(gco, 'linestyle'), '-') % ����ѡ���߲Ž��в���������*
    
    % ����ߴ���ѡ��״̬����λ�ڵ�ǰcell�ڣ���ȡ��ѡ��
    if strcmp(get( gco, 'Selected' ), 'on') && ~isempty(GT_move{t}{global_x,global_y}) && GT_move{t}{global_x,global_y}(end)==str2double(get(gco,'DisplayName'))
        
        set( gco, 'Selected', 'off' );
        GT_move{t}{global_x,global_y} = GT_move{t}{global_x,global_y}(1:end-1); % ɾȥcsp�е����һ��
        if isempty( GT_move{t}{global_x,global_y} )
            msgbox('CSP�����');
        end
        
    else % ���δѡ�У���ѡ�в���¼
        set( gco, 'Selected', 'on' );
        % д�뵱ǰλ��
        if all(GT_move{t}{global_x,global_y}~=str2double(get(gco, 'DisplayName')))
            GT_move{t}{global_x,global_y} = [ GT_move{t}{global_x,global_y}, str2double(get(gco, 'DisplayName')) ];
        end
        
    end
    
end
   





end

