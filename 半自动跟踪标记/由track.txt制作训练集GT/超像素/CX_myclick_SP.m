function CX_myclick_SP( hObject, eventdata, handles )

global tmp_label2e global_x global_y last_click;

if ~strcmp( get( gco, 'Type' ), 'line') % ����ѡ���߲Ž��в��������Ǳ���
%     msgbox('ûѡ�ж�������ѡ');
else
    % 2015.8.27 �¸Ķ���ͨ��ʶ��ǰѡ�ж��󣨵���ߣ��������б�ǣ���˿�����ȷ���merge�¼�
    % 2015.10.9 �¸Ķ���ͨ��ȫ�ֱ��� last_click ����ֹ��ѡ�ߺ�ѡ��
    % dot_or_lineΪ1ʱ˵���յ����һ����
    isdot = isequal( get( gco, 'marker' ), '*');
    
    % ��겻��1��ȴѡ�е㣬������1��ȴѡ���ߣ������������
    if global_y~=1 && isdot && ~strcmp(get(gco, 'DisplayName'), last_click) 
        msgbox('ȷ������밴�ո����ѡ����һ����!');
        return
    end
    if global_y==1 && ~isdot
        msgbox('����ѡ���!');
        return
    end
        
    if strcmp(get(gco, 'DisplayName'), last_click) && strcmp(get(gco, 'Selected'),'on')
        % ѡ������һ����ѡ�ģ������ɾ��
        set( gco, 'Selected', 'off' );
        global_y = global_y - 1; % ����˺�һ�񣬽�������ѡȡ
        if global_y==1
            msgbox('ɾ��һ�����ӣ�');
        end    
    else
        % ����ѡ����һ��ѡ�ģ�����м�¼
        set( gco, 'Selected', 'on' );
        tmp_label2e(global_x, global_y) = str2double( get(gco, 'DisplayName') ); % ��¼
        % ѡ�����������𣬼���ѡ����һ��
        global_y = global_y + 1;
        
        
        if global_y>size(tmp_label2e,2)
            msgbox('ѡ���Superpixel��Ŀ�ѵ����ޣ��벻Ҫ��ѡ��');
        end
    end
    
    last_click = get(gco, 'DisplayName'); % ����last_click

end
    
