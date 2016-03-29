function [PRE, REC, FM, Tcount, Pcount] = cal_move_prec_rec_Fm( man_track_txt, res_track_txt )


%% 读入tif图片和track。txt
load(man_track_txt);
load(res_track_txt);
man_dir_path = man_track_txt(1:end-13);
man_dir = dir([man_dir_path,'*.tif']);
res_dir_path = res_track_txt(1:end-13);
res_dir = dir([res_dir_path,'*.tif']);

% 进行数据与处理
% truth 数目 predict数目
man_track(:,2:3) = man_track(:,2:3) + 1;
res_track(:,2:3) = res_track(:,2:3) + 1;
Tcount = sum(man_track(:,3) - man_track(:,2));
Pcount = sum(res_track(:,3) - res_track(:,2));

%% 统计TP数目
TP = 0;
frame = max(man_track(:,3));
for tt=1:frame-1
    % 读入4张图片
%     fprintf('正在处理 %d - %d 帧...\n',tt,tt+1);
    gt1 = imread([man_dir_path, man_dir(tt).name]);
    res1 = imread([res_dir_path, res_dir(tt).name]);
    gt2 = imread([man_dir_path, man_dir(tt+1).name]);
    res2 = imread([res_dir_path, res_dir(tt+1).name]);
    
    % 处理前2帧
    Lgt1 = regionprops(gt1);
    Lres1 = regionprops(res1);
    g1 = cat(1, Lgt1.Centroid);
    r1 = cat(1, Lres1.Centroid);
    distance = dist2(g1,r1); % 把gt匹配到最近的res标签上
    matches1 = zeros(numel(Lgt1),1);
    for i=1:numel(Lgt1)
        if Lgt1(i).Area==0
            matches1(i) = 0;
            continue
        end
        [~,matches1(i)] = min(distance(i,:));
        distance(:,matches1(i)) = inf; % 将用过的标签的距离置为inf，防止别的再继续匹配上
    end
    % 处理后2帧
    Lgt2 = regionprops(gt2);
    Lres2 = regionprops(res2);
    g2 = cat(1, Lgt2.Centroid);
    r2 = cat(1, Lres2.Centroid);
    distance = dist2(g2,r2); % 把gt匹配到最近的res标签上
    matches2 = zeros(numel(Lgt2),1);
    for i=1:numel(Lgt2)
        if Lgt2(i).Area==0
            matches2(i) = 0;
            continue
        end
        [~,matches2(i)] = min(distance(i,:));
        distance(:,matches2(i)) = inf;
    end
    % 统计tp的数目
    for j=1:numel(matches1)
        if j>numel(matches2) % 1中多于2，则不是迁移
            break
        end
        if matches1(j)~=0 && matches1(j)==matches2(j)
            TP = TP + 1;
        end
    end
end

%% 最终结果
fprintf('\n迁移3值如下：\n');
PRE = TP/Pcount; fprintf('\nPrecision:\t%f\n',PRE);
REC = TP/Tcount; fprintf('Recall:\t\t%f\n',REC);
FM = 2/(1/PRE+1/REC); fprintf('F-measure:\t%f\n',FM);
fprintf('Truth number:\t%d\n', Tcount);
fprintf('Predict number:\t%d\n', Pcount);










