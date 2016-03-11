%######################################
%
% 2015.6.6 CX on desk
% ���ã�����������ڴ�GT TRA�еõ��Ҷȱ�ǩ���ҵ���Բ��ŵĶ�Ӧ��ϵ�����ð��Զ���� 
% ���ݴ洢�����õ��Ķ�Ӧ���󱣴�Ϊ Label_to_Ellipse.mat 
% ������ϵ������ CX_myclick ��ΪGUI��������Ӧ����������ʵ��ѡ��/ȡ��
%
%######################################

clear;close all;
[ ~, trackpath ] = getpath( 'training' );

last = max(strfind(trackpath, '\'));
gtpath = [trackpath(1:last+2), '_GT\TRA\'];

gt_dir = dir([gtpath, '*.tif']);
fig_dir = dir([trackpath,'\�����ͼ\*.fig']);
output_dir = dir([trackpath,'\GT\label_and_e\*.fig']);
center_gt_path = [trackpath, '\GT\center_gt.mat'];

frame = numel(fig_dir);

%% ���� TRA �б�ǩ*����Բ���λ��ͼ

if exist(center_gt_path, 'file')
    load(center_gt_path);
else
    center_gt = cell(frame,1);
    stats = cell(frame,1);
    % �ȼ���label�����ݣ�������������Ͱ볤��
    tic
    for t=1:frame
        gt = imread([ gtpath, gt_dir(t).name ]);
        stats{t} = regionprops(gt, 'Centroid', 'MajorAxisLength'); % ������Բ��һЩ��Ϣ
        center_gt{t} = cat(1,stats{t}.Centroid);
    end
    toc
    save(center_gt_path, 'stats','center_gt');
end

% ---------------------------------------- %
if 0 % �Ƿ����label2eͼ��
    plot_label2e( stats, fig_dir );
end

%% ����ԭʼ��Բ��Ϣ
load([ trackpath, '\Pair\Pre_data_new.mat'], 'Ellipse','n');

% ���� label �� ellipse �Ķ�Ӧ��ϵ����
if ~exist([ trackpath, '\GT\Label_to_Ellipse.mat'],'file')
    label2e = cell(frame,1); % label2e ����Ϊlabel��ellipse�Ķ�Ӧ��ϵ
else
    load([ trackpath, '\GT\Label_to_Ellipse.mat']);
end

center_e = cell(frame,1); % ������Բ�����ĵ�
distance = cell(frame,1); % ����������Բ��*�ľ���
% ���ӱ��ʱ����Ҫ�Ծ����������
screen_size = get(0,'ScreenSize');

%% ���а��Զ����

for t=92:frame % ��frame���У���Ŀǰ�����û�еĽ��б��
    for j=1:n(t)
        center_e{t}(j,1) = Ellipse{t}{j}.x0;
        center_e{t}(j,2) = Ellipse{t}{j}.y0;
    end
    
    %% ����˵ǰ�������Զ���ǣ�ÿ��*�ҵ����Լ��������Բ
    distance{t} = dist(center_e{t}, center_gt{t}'); 
    
    for label=1:size(center_gt{t},1) % labelΪGT��ǩ�Ҷ�ֵ
        if isnan( distance{t}(1,label) ) % �������ΪNaN�����ӦΪNaN
            label2e{t}(label,1) = NaN;
            continue;   
        end
        % �ҳ�����*���������Բ���,��Ϊ�Ǿ������������һ��*Ψһ��Ӧһ����Բ����һ����Բ���ܶ�Ӧ2��*
        j = find( distance{t}(:,label) == min(distance{t}(:,label)) ); 
        % ���*������Բ�ڸ�����������ԲΪ����˵ǰ������ֱ�Ӷ�Ӧ��
        if distance{t}(j,label)<= Ellipse{t}{j}.b + stats{t}(label).MajorAxisLength && Ellipse{t}{j}.num_hypoth ==1
            label2e{t}(label,1) = j;
        else
            label2e{t}(label,1) = 0; % ������Ϊ0��˵������Ҫ�˹����
        end
    end
    
    %% �ҳ�δ��label��Ӧ�ϵ���Բ�������龰�Ͷ��˵ǰ���������ֶ����
    e_not_labed = setdiff( 1:numel(Ellipse{t}), label2e{t}(:,1) );
    label_not_attached = find(label2e{t}==0);
    if isempty(e_not_labed) && isempty(label_not_attached)
        continue;
    end
    
    % �򿪵�ǰͼƬ����������Բ��ǳɺ�ɫ
    disp(['��ǰ���ڴ��� ',fig_dir(t).name]);
    openfig( [ trackpath, '\GT\label_and_e\', fig_dir(t).name] );
    set(gcf, 'Position', screen_size);
    % ȫ����ʾ
%     screen_size = get(0,'ScreenSize');
%     set(gcf,'Position',screen_size);
    h = get(gca, 'children');
%     h = findobj(h, 'Type','Line');
    % ���ð�������
    global tmp_label2e global_x global_y last_click;
    tmp_label2e = zeros(50,2); % ��������ڴ��label��e�ı��
    global_x = 1;global_y = 1;last_click = 0;
    % ����Ϊ����->��ʾ
    set(h, 'ButtonDownFcn', @CX_myclick); 
    % ��Ҫ�ҳ�������Բ�ľ�������ܽ��б��
    % ------------------------------------------------------------------- %
    % ��û���Զ�����ϵ���Բ����Ϊ��ɫ��*�����Ϊ��ɫ���Լ�ǿ�Ӿ�����
    % 1��������ɫ��
    for ind = 1:numel(label_not_attached)
        h_label = findobj(h, 'color', 'w', 'DisplayName', num2str(label_not_attached(ind)) );
        set(h_label, 'color', 'b');
    end
    % 2�����ƺ�ɫ��Բ
    for ind_j=1:numel(e_not_labed)
        h_e = findobj(h, 'color', 'g', 'DisplayName', num2str(e_not_labed(ind_j)) );
        set(h_e, 'color', 'r'); % ��ɫ
    end
    % ------------------------------------------------------------------- %
    
    % ������ͣ���ȴ�������ͼ��
    keyin = input('  �ڿ���̨�ϰ�enter��������һ֡\n', 's');
    if strcmp(keyin, '')
        close(gcf);
    end

    % �����Ǻ�� tmp_label2e ����ȥ��0��nan
    tmp_label2e = tmp_label2e(tmp_label2e(:,1)~=0, :);
    tmp_label2e = tmp_label2e(~isnan(tmp_label2e(:,1)), :);
    % ����Ӧ��ϵ���� label2e ����
    label2e{t}(tmp_label2e(:,1),1) = tmp_label2e(:,2);
    clear global;
    
    if 1
        disp('�������� label2e.mat ��');
        save([ trackpath, '\GT\Label_to_Ellipse.mat'], 'label2e');
    end
    
end
%
% �˹����ʹ�÷����� ��ɫ����Բ��û����Ӧ�ϵ���Բ���ȵ���׵㣬�ٵ����Բ���Ϳ��Խ����Ӧ��
%                   �������ˣ������ٵ��һ�½���ȡ��
%


    
    
    





















        


