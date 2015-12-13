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


% ����ͼƬ
img = imread(raw_pic);
bw = imread(bw_pic);
bw = imfill(bw, 'holes');

img(~bw) = 0; % bwǷ�ָ�Ϊ0�ľ�ǿ��ʹԭʼͼƬҲΪ0
imga = imadjust(img); % ��ת��double��ʽ imshow(imga)
% ע�⣺������˳�����Ҫ��������=0��adjustЧ�����ã���Ϊǰ���������������ڷָ�

imgs = im2uint8(im2double(imga)); % ֻ֧��uint8��Ϊ���� imshow(imgs)
[labels, numlabels] = slicomex(imgs, nsp); % numlabels is the same as number of superpixels

if isa(imga, 'uint16')
    col = 65535; % ��ɫҪע������ͼƬ��ά��
elseif isa(imga, 'uint8')
    col = 255;
elseif isa(imga, 'double')
    col = 1;
end

% ����ΪͼƬ��ʾЧ����չʾ�ָ��ߣ�
maskim = drawregionboundaries(labels, imga, col);
% imshow(maskim);

%% ʹ��mask����superpixel��Ч��
[h,w] = size(labels);
raw_stats = regionprops(labels, img, 'PixelValues');
for row=1:numel(raw_stats)
    % ����ɫ�����ؿ��label��Ϊ0��ֻ�ж����ĵ��Ƿ�Ϊ��ɫ����������2015.12.8
%     xyco = round([stats(row).Centroid(1), stats(row).Centroid(2)]);
%     if bw(xyco(2),xyco(1))==0 % ע������Ҫ��д
%         labels(labels==row) = 0;
%     end
    
    % ���ǲ���ͳ�ƺ�ɫ������ķ����ȽϺã�2015.12.10��Ϊͳ�ƷǺڵ�����ķ�����
%     percent = sum(raw_stats(row).PixelValues==0)/numel(raw_stats(row).PixelValues);
%     if percent>=0.8 % ��ɫ��������80%��������0
%         labels(labels==row) = 0;
%     end
    if sum(raw_stats(row).PixelValues~=0)<rm_small % �Ǻڵ�̫�٣���ȫ����0
        labels(labels==row) = 0;
    end
end

%% ��mask�к�ɫ���ָ�ֵ��labels������ᵼ��ԭ��1��label����Ϊ2�������Ҫ���������label��ֵ
% ������������Ч����Ƿ�ָ����Ϣ
labels(~bw) = 0;

new_l = numlabels-1; % ԭ����label��Ŀ����1����Ϊ������0���µı�Ŵ����￪ʼ

has_small = true; % ��ʶ���������ж��Ƿ���С����
while has_small % ���У��򲻶�ִ��ɾ��������ɾ�����ܻᵼ���������仯�����Ҫ��2����ͬ��ѭ����
    has_small = false;
    
    for row=unique(labels)'
        if row==0 % ����������
            continue;
        end
        % ����Ƿ�Ϊ������
        bb = zeros(size(labels));
        bb(labels==row) = 1; % imshow(bb) �����鿴��label����
        [L, n_region] = bwlabel(bb, 4);
        if n_region==1 % �������һ����ͨ����������
            continue
        end
        % �Զ����������������и��±�ǩ
        for j=2:n_region
            new_l = new_l + 1;
            labels(L==j) = new_l;
        end
    end
    
    for row=unique(labels)'
        % ---------------------------------------- %
        % ȥ��̫С��label
        if sum(sum(labels==row))<rm_small 
            has_small = true;
            labels(labels==row) = 0;
        end
        % ---------------------------------------- %
    end

end
disp(['  ͨ��mask��label��Ŀ��',num2str(numlabels-1),'���ӵ���',num2str(new_l)]);
    


% ʹ��maskʹ�ô���label��Ϊ0��Ҫ��0����
new_labels = zeros(size(labels));
seql = unique(labels);
for n=1:numel(seql)
    % ��ԭ���ı�ǩ��0��ʼ��˳������
    new_labels(labels==seql(n)) = n-1;
end

% imagesc(new_labels)
% axis image
% ��ʾlabelͼ��
% figure
RGB_label = label2rgb(new_labels,@hsv,'k','shuffle');
% imshow(RGB_label);

%% �ҳ����ڵĳ�����
label_nearby = cell(numel(seql)-1,1); % seql��һ����0�������ڱ�ǩ
for row=1:h-1
    for j=1:w-1
        this = new_labels(row,j);
        r_neighbor = new_labels(row,j+1); % �ұߵ��ھ�
        b_neighbor = new_labels(row+1,j); % �±ߵ��ھ�
        if this==0 || r_neighbor==0 || b_neighbor==0
            continue;
        end
        % ����һ���ұ��ھ�
        if  r_neighbor~=this && ~any(label_nearby{this}==r_neighbor)
            label_nearby{this} = [ label_nearby{this}, r_neighbor ];
        end
        % ����һ���±��ھ�
        if b_neighbor~=this && ~any(label_nearby{this}==b_neighbor)
            label_nearby{this} = [ label_nearby{this}, b_neighbor ];
        end
    end
end
        
% ���label����֮������ڹ�ϵ
for row=1:numel(label_nearby)
    for j=label_nearby{row}
        if ~any(label_nearby{j}==row)
            label_nearby{j} = [label_nearby{j}, row];
        end
    end
end

%% ���ɾ����˵

n_fore = 0;
foreground = cell(1); % �к�Ϊǰ����ţ��޹��򣩣�����Ϊǰ���ڵ�label��
global cluster

for row=1:numel(label_nearby)
    n_fore = n_fore + 1;
    % �ж�֮ǰ��foreground�Ƿ��������һ��
    res = cellfun(@(x) any(x==row),foreground,'un',1);
    if any(res==1) % ���������򲻴�������
        n_fore = n_fore - 1;
        continue;
    end
    
    cluster = [];
    % �������� add_cluster
    add_into_cluster( label_nearby, row) 
    foreground{n_fore,1} = cluster;
end

%% superpixels �к�Ϊǰ����ţ��洢��˵��Ϣ������ellipse��
bsp_stats = regionprops(new_labels, 'pixellist');
superpixels =  cell(size(foreground));
for row=1:numel(superpixels) % i����ǰ�����
    % �������е�������
    fore = foreground{row};
    [ label_zuhe, flag_zuhe ] = generate_csp( fore );
    % �������Щ�ǲ����ڵģ�Ҫȥ��
    [ label_zuhe flag_zuhe ] = delete_uncombined_sp( fore, label_nearby, label_zuhe, flag_zuhe, bsp_stats, new_labels );
    % ��ÿ����˵������Ϣ���൱����Բ�е���Ϣ��
    superpixels = calculate_sp_info( row, label_zuhe, flag_zuhe, bsp_stats, superpixels );
end

% �� superpixels ��ֱ��һ��
n_in_fore = cellfun(@numel, superpixels); % ÿ��ǰ���еļ�˵��
n_all = sum(n_in_fore); % ����ͼ�ļ�˵��
SP = cell(n_all,1);
for row=1:numel(superpixels)
    i_s = sum(n_in_fore(1:row-1))+1; % �������
    i_e = sum(n_in_fore(1:row)); % �����յ�
    
    SP(i_s:i_e) = superpixels{row};
end
        












