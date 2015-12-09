%######################################
%
% 2015.12.9 CX on desk
% ���ã�����������ڴ�GT TRA�еõ��Ҷȱ�ǩ���ҵ���Բ��ŵĶ�Ӧ��ϵ�����ð��Զ���� 
% ���ݴ洢�����õ��Ķ�Ӧ���󱣴�Ϊ Label_to_Ellipse.mat 
% ������ϵ������ CX_myclick ��ΪGUI��������Ӧ����������ʵ��ѡ��/ȡ��
%
%######################################

clear;close all;
[ segpath, trackpath ] = getpath( 'training' );

last = max(strfind(segpath, '\'));
gtpath = [segpath(1:last+2), '_GT\TRA\'];
gt_dir = dir([gtpath, '*.tif']); % gtͼƬ��λ��

frame = numel(gt_dir);
frame = 10;

% fig������ǵ�ͼƬ
fig_path = [trackpath,'\GT\label_and_e\'];
mkdir(fig_path)
fig_dir = dir([fig_path, '\*.fig']);

%% ���� TRA �б�ǩ*����Բ���λ��ͼ
center_gt_path = [trackpath, '\GT\center_gt.mat'];
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
    plot_label2e_SP( stats );
end

%% ����ԭʼ��Բ��Ϣ
load([ trackpath, '\Pair\Pre_data_new.mat'], 'SuperPixel','n');

% ���� label �� ellipse �Ķ�Ӧ��ϵ����
if ~exist([ trackpath, '\GT\Label_to_Ellipse.mat'],'file')
    label2e = cell(frame,1); % label2e ����Ϊlabel��ellipse�Ķ�Ӧ��ϵ
else
    load([ trackpath, '\GT\Label_to_Ellipse.mat']);
end

center_sp = cell(frame,1); % ����SP�����ĵ�
distance = cell(frame,1); % ��������SP��*�ľ���
% ���ӱ��ʱ����Ҫ�Ծ����������
screen_size = get(0,'ScreenSize');

%% ���а��Զ����
global tmp_label2e global_x global_y last_click;

for t=3:frame % ��frame���У���Ŀǰ�����û�еĽ��б��
    for j=1:n(t)
        center_sp{t}(j,:) = SuperPixel{t}{j}.centroid;
    end
    %% ����˵ǰ�������Զ���ǣ�ÿ��*�ҵ����Լ��������Բ
    distance{t} = dist(center_sp{t}, center_gt{t}'); 
    
    for label=1:size(center_gt{t},1) % labelΪGT��ǩ�Ҷ�ֵ
        if isnan( distance{t}(1,label) ) % �������ΪNaN�����ӦΪNaN
            label2e{t}(label,1) = NaN;
            continue;   
        end
        % �ҳ�����*���������Բ���,��Ϊ�Ǿ������������һ��*Ψһ��Ӧһ����Բ����һ����Բ���ܶ�Ӧ2��*
        j = find( distance{t}(:,label) == min(distance{t}(:,label)) ); 
        % ���*������Բ�ڸ�����������ԲΪ����˵ǰ������ֱ�Ӷ�Ӧ��
        if distance{t}(j,label)<= stats{t}(label).MajorAxisLength && SuperPixel{t}{j}.num_hypoth ==1
            label2e{t}(label,1) = j; % ��һ�д��SP���
            label2e{t}(label,2) = SuperPixel{t}{j}.label; % �ڶ��д��SP��ǩ
        else
            label2e{t}(label,1) = 0; % ������Ϊ0��˵������Ҫ�˹����
        end
    end
    
    %% �ҳ�δ��label��Ӧ�ϵ���Բ�������龰�Ͷ��˵ǰ���������ֶ����
    
    disp(['��ǰ���ڴ��� ',fig_dir(t).name]);
    openfig( [ trackpath, '\GT\label_and_e\', fig_dir(t).name] );
    % ȫ����ʾ
    set(gcf, 'Position', screen_size);
    
    h = get(gca, 'children');
    % ���ð�������
    
    tmp_label2e = zeros(50,6); % ��������ڴ��label��e�ı��
    global_x = 1;
    global_y = 1;
    last_click = '';
    % ����Ϊ����->��ʾ
    set(h, 'ButtonDownFcn', @CX_myclick_SP); 
    set(gcf, 'keypressfcn', @pressBlank_SP); % ���ո��������')
    % ��Ҫ�ҳ�������Բ�ľ�������ܽ��б��
    
    % ------------------------------------------------------------------- %
    % ��û���Զ�����ϵ���Բ����Ϊ��ɫ��*�����Ϊ��ɫ���Լ�ǿ�Ӿ�����
    % 1��������ɫ��
    label_not_attached = find(label2e{t,1}==0);
    for lab=label_not_attached'
        h_label = findobj(h, 'marker','*','DisplayName',num2str(lab) );
        set(h_label, 'color', 'b');
    end
    % 2�����ƺ�ɫԲȦ
    bsp = cell2mat( cellfun(@(x) x.label, SuperPixel{t},'un',0)' );
    maxbsp = max(bsp); % �ҵ����ͼ������basic sp��Ҳ����basic sp������Ŀ��
    sp_not_labed = setdiff( 1:maxbsp, label2e{t}(:,2) );
    for jj=sp_not_labed
        h_e = findobj(h, 'color','g','DisplayName',num2str(jj) );
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
%     label2e{t}(:,2) = []; % ���label2e�ĵڶ��У�ֻ�ڻ���ɫԲʱ���ã�
    for i_row=1:size(tmp_label2e,1)
        row = tmp_label2e(i_row,:);
        tmp = row(row~=0);
        this_sp = tmp(2:end);
        flag = cellfun(@(x) isequal(sort(this_sp), sort(x.label)), SuperPixel{t});
        assert(any(flag==1)) % ����һ�����ҵ��Ǹ�label���
        label2e{t}(row(1),1) = find(flag);
        label2e{t}(row(1),2:numel(this_sp)+1) = this_sp;
    end
    
    if 1
        disp('�������� label2e.mat ��');
        save([ trackpath, '\GT\Label_to_Ellipse.mat'], 'label2e');
    end
    
end
%
% �˹����ʹ�÷����� ��ɫ����Բ��û����Ӧ�ϵ���Բ���ȵ���׵㣬�ٵ����Բ���Ϳ��Խ����Ӧ��
%                   �������ˣ������ٵ��һ�½���ȡ��
%


    
    
    





















        


