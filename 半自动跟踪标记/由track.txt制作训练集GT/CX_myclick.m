function CX_myclick( hObject, eventdata, handles )

global tmp_label2e global_x global_y last_click;

if ~strcmp( get( gco, 'Type' ), 'line') % ����ѡ���߲Ž��в��������Ǳ���
%     msgbox('ûѡ�ж�������ѡ');
else
    % 2015.8.27 �¸Ķ���ͨ��ʶ��ǰѡ�ж��󣨵���ߣ��������б�ǣ���˿�����ȷ���merge�¼�
    % 2015.10.9 �¸Ķ���ͨ��ȫ�ֱ��� last_click ����ֹ��ѡ�ߺ�ѡ��
    % dot_or_lineΪ1ʱ˵���յ����һ����
    isdot = isequal( get( gco, 'marker' ), '*'); 
    if (global_y==1 && isdot) || (global_y==2 && ~isdot) % �����1����ѡ�е���߹����2����ѡ���ߣ����¼
        set( gco, 'Selected', 'on' );
        last_click = get(gco, 'DisplayName'); % �յ�����Ŀ��
        tmp_label2e(global_x,global_y) = str2double( last_click );
        % ѡ�����������𣬼���ѡ����һ��
        global_y=global_y+1;
        if global_y==3 % ˵��������¼�����
            msgbox('���һ�����ӣ�'); 
            global_y=1;
            global_x=global_x+1;
        end
    else % ������������1ѡ����/���2ѡ�е㣬������ɾ������
        if ~strcmp(get(gco, 'DisplayName'), last_click) 
            % ���ǵ��֮ǰ��ѡ�еģ��򲻲���
        else
            set( gco, 'Selected', 'off' );
                if global_y==2 % ����˺�һ�񣬽�������ѡȡ
                    global_y=1;
                    msgbox('ɾ��һ�����ӣ�');
                else
                    global_y=2;
                    global_x=global_x-1;   
                    last_click = num2str(tmp_label2e(global_x,1)); % ��ѡ�е�ҲҪ����һ��
                end
                if global_x>0 % ��һ�ξ�ѡ��ʱglobal_x�����0����Ҫ�����ж�
                    tmp_label2e(global_x,global_y) = 0; % �˺�һ�����û�н�������ѡȡ����λ���ϵ������ڴ��ڣ����Ҫ������ 0
                else
                    error('������ѡ�����ѡ����Բ��');
                end
        end
    end
    
% ��������Ϊ�ɰ棬�޷����merge
%     if strcmp(get( gco, 'Selected' ), 'off') % ����ߴ���δѡ��״̬����ѡ����
%         set( gco, 'Selected', 'on' );
%         tmp_label2e(global_x,global_y) = str2double( get(gco, 'DisplayName') );
%         % ѡ�����������𣬼���ѡ����һ��
%         global_y=global_y+1;
%         if global_y==3
%             global_y=1;
%             global_x=global_x+1;
%         end
%     else % ����ѱ�ѡ�У���ȡ��ѡ��״̬         
%         set( gco, 'Selected', 'off' );
%         if global_y==2 % ����˺�һ�񣬽�������ѡȡ
%             global_y=1;
%         else
%             global_y=2;
%             global_x=global_x-1;
%         end
%         tmp_label2e(global_x,global_y) = 0; % �˺�һ�����û�н�������ѡȡ����λ���ϵ������ڴ��ڣ����Ҫ������ 0
%     end

end
    
 
end