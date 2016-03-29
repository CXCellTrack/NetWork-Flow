% 统计分割精度

clear
% 载入标准答案的中心点数据
load('E:\datasets\first_edition\training_datasets\N2DH-SIM\05_0-00_track\GT\center_gt.mat');
load('E:\datasets\first_edition\training_datasets\N2DH-SIM\05_0-00_track\Pair\Pre_data_new.mat','n');
contour_path = 'E:\datasets\first_edition\training_datasets\N2DH-SIM\05_0-00_seg\FOI提取轮廓\';

frame = numel(n);
contour_dir = dir([ contour_path, '*.tif' ]);

%% 进行一一对应，对应不上的可计算prec和recall
center_e = cell(frame,1);
distance = cell(frame,1);
TP = zeros(frame,1);   % 前景和细胞对应个数
Tnum = zeros(frame,1); % 实际细胞个数
Pnum = zeros(frame,1); % 分割得到的前景个数

for t=1:frame % 对frame中有，但目前标记中没有的进行标记
    disp(['处理第',num2str(t),'帧...']);
    center_gt{t} = center_gt{t}(~isnan(center_gt{t}(:,1)),:); % 去掉center中为nan的那些行
    
    im = imread([contour_path, contour_dir(t).name]); % imshow(im)
    im = bwareaopen(im, 50);
    stats = regionprops(im, 'centroid','MajorAxisLength');
    center_region = cat(1, stats.Centroid);
    % 为每个前景找到距离最近的*点，那么每个前景至多有一个*点，此时为TP；若找不到*点，为FP；多出来的*点为FN
    distance{t} = dist(center_region, center_gt{t}'); 
    for re=1:numel(stats)
        [thisD, label] = min(distance{t}(re,:)); 
        if thisD<=stats(re).MajorAxisLength % 如果距离小于区域圆心（即*点位于区域内的快速近似判断）
            TP(t) = TP(t) + 1;
        end
    end
    
    [ Pnum(t), Tnum(t) ] = size(distance{t});
    
end
     
precision = sum(TP)/sum(Pnum);
recall = sum(TP)/sum(Tnum);
F_m = precision*recall*2/(precision+recall);
% 打印输出
fprintf('precision:\t%f\nrecall:\t\t%f\nF_m:\t\t%f\n',precision,recall,F_m);
    
    
    
    
    
    
    
    
    
    
    
    
    
    