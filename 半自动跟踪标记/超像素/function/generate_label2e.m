function generate_label2e()

[ segpath, trackpath ] = getpath( 'training' );

last = max(strfind(segpath, '\'));
gtpath = [segpath(1:last+2), '_GT\TRA\'];
gt_dir = dir([gtpath, '*.tif']); % gtͼƬ��λ��

frame = numel(gt_dir);

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

%% �Զ���Ǽ򵥵ĳ�����

for t=1:frame % ��frame���У���Ŀǰ�����û�еĽ��б��
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
    
end

if 1
    disp('�������� label2e.mat ��');
    save([ trackpath, '\GT\Label_to_Ellipse.mat'], 'label2e');
end
    
    
    





















        


