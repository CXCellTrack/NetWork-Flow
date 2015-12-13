%======================================================================
%SLICO demo
% Copyright (C) 2015 Ecole Polytechnique Federale de Lausanne
% File created by Radhakrishna Achanta
% Please also read the copyright notice in the file slicmex.c 
%======================================================================
%Input:
%[1] 8 bit images (color or grayscale)
%[2] Number of required superpixels (optional, default is 200)
%
%Ouputs are:
%[1] labels (in raster scan order)
%[2] number of labels in the image (same as the number of returned
%superpixels
%
%NOTES:
%[1] number of returned superpixels may be different from the input
%number of superpixels.
%[2] you must compile the C file using mex slicmex.c before using the code
%below
%----------------------------------------------------------------------
% How is SLICO different from SLIC?
%----------------------------------------------------------------------
% 1. SLICO does not need compactness factor as input. It is calculated
% automatically
% 2. The automatic value adapts to the content of the superpixel. So,
% SLICO is better suited for texture and non-texture regions
% 3. The advantages 1 and 2 come at the cost of slightly poor boundary
% adherences to regions.
% 4. This is also a very small computational overhead (but speed remains
% almost as fast as SLIC.
% 5. There is a small memory overhead too w.r.t. SLIC.
% 6. Overall, the advantages are likely to outweigh the small disadvantages
% for most applications of superpixels.
%======================================================================
function [ SP, new_labels, maskim, RGB_label ] = CX_Generate_SP_in_1( raw_pic, bw_pic, nsp, rm_small )


% 读入图片
img = imread(raw_pic);
bw = imread(bw_pic);
bw = imfill(bw, 'holes');

img(~bw) = 0; % bw欠分割为0的就强制使原始图片也为0
imga = imadjust(img); % 先转成double形式 imshow(imga)
% 注意：操作的顺序很重要，发现先=0再adjust效果更好，因为前景更加明亮，利于分割

imgs = im2uint8(im2double(imga)); % 只支持uint8作为输入 imshow(imgs)
[labels, numlabels] = slicomex(imgs, nsp); % numlabels is the same as number of superpixels

if isa(imga, 'uint16')
    col = 65535; % 颜色要注意输入图片的维数
elseif isa(imga, 'uint8')
    col = 255;
elseif isa(imga, 'double')
    col = 1;
end

% 仅作为图片演示效果（展示分割线）
maskim = drawregionboundaries(labels, imga, col);
% imshow(maskim);

%% 使用mask改善superpixel的效果
[h,w] = size(labels);
raw_stats = regionprops(labels, img, 'PixelValues');
for row=1:numel(raw_stats)
    % 将黑色超像素块的label变为0（只判断中心点是否为黑色）这样有误差！2015.12.8
%     xyco = round([stats(row).Centroid(1), stats(row).Centroid(2)]);
%     if bw(xyco(2),xyco(1))==0 % 注意坐标要反写
%         labels(labels==row) = 0;
%     end
    
    % 还是采用统计黑色点比例的方法比较好（2015.12.10改为统计非黑点个数的方法）
%     percent = sum(raw_stats(row).PixelValues==0)/numel(raw_stats(row).PixelValues);
%     if percent>=0.8 % 黑色比例大于80%，则将其置0
%         labels(labels==row) = 0;
%     end
    if sum(raw_stats(row).PixelValues~=0)<rm_small % 非黑点太少，则全部置0
        labels(labels==row) = 0;
    end
end

%% 将mask中黑色部分赋值给labels，但这会导致原本1个label被分为2个，因此要给多出来的label赋值
% 这样做可以有效利用欠分割的信息
labels(~bw) = 0;

new_l = numlabels-1; % 原本的label数目（减1是因为不包括0）新的编号从这里开始

has_small = true; % 标识符，用来判断是否有小区域
while has_small % 若有，则不断执行删除操作（删除可能会导致区域发生变化，因此要用2个不同的循环）
    has_small = false;
    
    for row=unique(labels)'
        if row==0 % 背景不处理
            continue;
        end
        % 检查是否为单区域
        bb = zeros(size(labels));
        bb(labels==row) = 1; % imshow(bb) 单独查看该label区域
        [L, n_region] = bwlabel(bb, 4);
        if n_region==1 % 如果还是一个联通区域，则跳过
            continue
        end
        % 对多区域的其他几块进行赋新标签
        for j=2:n_region
            new_l = new_l + 1;
            labels(L==j) = new_l;
        end
    end
    
    for row=unique(labels)'
        % ---------------------------------------- %
        % 去除太小的label
        if sum(sum(labels==row))<rm_small 
            has_small = true;
            labels(labels==row) = 0;
        end
        % ---------------------------------------- %
    end

end
disp(['  通过mask将label数目从',num2str(numlabels-1),'增加到了',num2str(new_l)]);
    


% 使用mask使得大量label变为0，要将0除掉
new_labels = zeros(size(labels));
seql = unique(labels);
for n=1:numel(seql)
    % 将原本的标签从0开始按顺序排列
    new_labels(labels==seql(n)) = n-1;
end

% imagesc(new_labels)
% axis image
% 显示label图像
% figure
RGB_label = label2rgb(new_labels,@hsv,'k','shuffle');
% imshow(RGB_label);

%% 找出相邻的超像素
label_nearby = cell(numel(seql)-1,1); % seql第一个是0，不属于标签
for row=1:h-1
    for j=1:w-1
        this = new_labels(row,j);
        r_neighbor = new_labels(row,j+1); % 右边的邻居
        b_neighbor = new_labels(row+1,j); % 下边的邻居
        if this==0 || r_neighbor==0 || b_neighbor==0
            continue;
        end
        % 增加一个右边邻居
        if  r_neighbor~=this && ~any(label_nearby{this}==r_neighbor)
            label_nearby{this} = [ label_nearby{this}, r_neighbor ];
        end
        % 增加一个下边邻居
        if b_neighbor~=this && ~any(label_nearby{this}==b_neighbor)
            label_nearby{this} = [ label_nearby{this}, b_neighbor ];
        end
    end
end
        
% 添加label互相之间的相邻关系
for row=1:numel(label_nearby)
    for j=label_nearby{row}
        if ~any(label_nearby{j}==row)
            label_nearby{j} = [label_nearby{j}, row];
        end
    end
end

%% 生成聚类假说

n_fore = 0;
foreground = cell(1); % 行号为前景编号（无规则），内容为前景内的label号
global cluster

for row=1:numel(label_nearby)
    n_fore = n_fore + 1;
    % 判断之前的foreground是否包含了这一行
    res = cellfun(@(x) any(x==row),foreground,'un',1);
    if any(res==1) % 若包含，则不处理这行
        n_fore = n_fore - 1;
        continue;
    end
    
    cluster = [];
    % 迭代函数 add_cluster
    add_into_cluster( label_nearby, row) 
    foreground{n_fore,1} = cluster;
end

%% superpixels 行号为前景编号，存储假说信息（类似ellipse）
bsp_stats = regionprops(new_labels, 'pixellist');
superpixels =  cell(size(foreground));
for row=1:numel(superpixels) % i代表前景编号
    % 生成所有的组合情况
    fore = foreground{row};
    [ label_zuhe, flag_zuhe ] = generate_csp( fore );
    % 组合中有些是不相邻的，要去掉
    [ label_zuhe flag_zuhe ] = delete_uncombined_sp( fore, label_nearby, label_zuhe, flag_zuhe, bsp_stats, new_labels );
    % 对每个假说计算信息（相当于椭圆中的信息）
    superpixels = calculate_sp_info( row, label_zuhe, flag_zuhe, bsp_stats, superpixels );
end

% 将 superpixels 拉直成一列
n_in_fore = cellfun(@numel, superpixels); % 每个前景中的假说数
n_all = sum(n_in_fore); % 整张图的假说数
SP = cell(n_all,1);
for row=1:numel(superpixels)
    i_s = sum(n_in_fore(1:row-1))+1; % 区域起点
    i_e = sum(n_in_fore(1:row)); % 区域终点
    
    SP(i_s:i_e) = superpixels{row};
end
        












