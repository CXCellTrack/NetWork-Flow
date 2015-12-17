function generate_label2e()

[ segpath, trackpath ] = getpath( 'training' );

last = max(strfind(segpath, '\'));
gtpath = [segpath(1:last+2), '_GT\TRA\'];
gt_dir = dir([gtpath, '*.tif']); % gt图片的位置

frame = numel(gt_dir);

% fig操作标记的图片
fig_path = [trackpath,'\GT\label_and_e\'];
mkdir(fig_path)
fig_dir = dir([fig_path, '\*.fig']);

%% 绘制 TRA 中标签*与椭圆间的位置图
center_gt_path = [trackpath, '\GT\center_gt.mat'];
if exist(center_gt_path, 'file')
    load(center_gt_path);
else
    center_gt = cell(frame,1);
    stats = cell(frame,1);
    % 先计算label的数据，包括中心坐标和半长轴
    tic
    for t=1:frame
        gt = imread([ gtpath, gt_dir(t).name ]);
        stats{t} = regionprops(gt, 'Centroid', 'MajorAxisLength'); % 计算椭圆的一些信息
        center_gt{t} = cat(1,stats{t}.Centroid);
    end
    toc
    save(center_gt_path, 'stats','center_gt');
end

% ---------------------------------------- %
if 0 % 是否绘制label2e图像
    plot_label2e_SP( stats );
end

%% 载入原始椭圆信息
load([ trackpath, '\Pair\Pre_data_new.mat'], 'SuperPixel','n');

% 载入 label 和 ellipse 的对应关系矩阵
if ~exist([ trackpath, '\GT\Label_to_Ellipse.mat'],'file')
    label2e = cell(frame,1); % label2e 矩阵为label和ellipse的对应关系
else
    load([ trackpath, '\GT\Label_to_Ellipse.mat']);
end

center_sp = cell(frame,1); % 保存SP的中心点
distance = cell(frame,1); % 用来保存SP到*的距离

%% 自动标记简单的超像素

for t=1:frame % 对frame中有，但目前标记中没有的进行标记
    for j=1:n(t)
        center_sp{t}(j,:) = SuperPixel{t}{j}.centroid;
    end
    
    %% 单假说前景进行自动标记，每个*找到离自己最近的椭圆
    distance{t} = dist(center_sp{t}, center_gt{t}'); 
    
    for label=1:size(center_gt{t},1) % label为GT标签灰度值
        if isnan( distance{t}(1,label) ) % 如果该列为NaN，则对应为NaN
            label2e{t}(label,1) = NaN;
            continue;   
        end
        % 找出距离*点最近的椭圆编号,因为是距离最近，所以一个*唯一对应一个椭圆；但一个椭圆可能对应2个*
        j = find( distance{t}(:,label) == min(distance{t}(:,label)) ); 
        % 如果*点在椭圆内附近，并且椭圆为单假说前景，则直接对应上
        if distance{t}(j,label)<= stats{t}(label).MajorAxisLength && SuperPixel{t}{j}.num_hypoth ==1
            label2e{t}(label,1) = j; % 第一列存放SP编号
            label2e{t}(label,2) = SuperPixel{t}{j}.label; % 第二列存放SP标签
        else
            label2e{t}(label,1) = 0; % 这个标记为0则说明其需要人工标记
        end
    end
    
end

if 1
    disp('保存结果到 label2e.mat 中');
    save([ trackpath, '\GT\Label_to_Ellipse.mat'], 'label2e');
end
    
    
    





















        


