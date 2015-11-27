

clear;close all

kthpath = 'E:\datasets\second_edition\competition_datasets\Fluo-N2DH-SIM+\01_KTH_RES\';
kthdir = dir([kthpath,'*.tif']);
savepath = 'E:\datasets\second_edition\competition_datasets\Fluo-N2DH-SIM+\01_KTH_RES_TRACK\';

%% 读入track.txt
KTH_track = load([ kthpath, 'res_track.txt']); % 载入KTH的track.txt
KTH_track(:,2:3) = KTH_track(:,2:3) + 1;

% 统计事件
divide = KTH_track(KTH_track(:,4)~=0,:);
other = KTH_track(KTH_track(:,4)==0,:);

%% 选择颜色
color = colormap(hsv);close('1'); 
color_numbel = numel(color)/3;
rng(1);
color = color(randperm(color_numbel),:);
nc = size(color,1);
inten2color = zeros(500,1); % 灰度映射到颜色索引上


for n=1:numel(kthdir)
    
    disp(['  处理图片 ',num2str(n),'...']);
    im = imread([kthpath, kthdir(n).name]);
    imd1 = zeros(size(im));imd2 = imd1;imd3 = imd1;
    % 统计im中的区域  
    stats = regionprops(im, 'Centroid');
    % 第一帧分配颜色
    if n==1 && numel(stats)<=nc
    	inten2color(1:numel(stats)) = 1:numel(stats); 
    end
    
    %% 进行上色
    for j=1:numel(stats)
        if inten2color(j)==0 % 说明是新出现的
            if j<=nc % 没用完表则继续对应
                inten2color(j) = j;
            else % 表用完了则随机产生一个
                inten2color(j) = randi(nc);
            end
        end
        imd1(im==j) = color(inten2color(j),1);
        imd2(im==j) = color(inten2color(j),2);
        imd3(im==j) = color(inten2color(j),3); 
        imd = cat(3, imd1,imd2,imd3);
    end
    imshow(imd); hold on;

    %% 统计分裂并标记
    ind = find(divide(:,2)-1==n);
    if ~isempty(ind)
        assert(mod(numel(ind),2)==0); % 必须为偶数个
        for j=ind'
            son = divide(j,1);
            father = divide(j,4);
%             inten2color(son) = inten2color(father); % 继承颜色(容易看不清楚)
            inten2color(son) = randi(nc); % 另外开一种颜色
            % 进行标记
            e.x0 = stats(father).Centroid(1);
            e.y0 = stats(father).Centroid(2);
            CX_plot_event( e, 'd');
        end
    end
    
    %% 统计出现并标记
    if n~=1
        ind = find(KTH_track(:,2)==n);
        if ~isempty(ind)
            for j=ind'
                tmpe = KTH_track(j,1);
                e.x0 = stats(tmpe).Centroid(1);
                e.y0 = stats(tmpe).Centroid(2);
                
                if any(divide(:,4)==KTH_track(j,4)) % 如果是子细胞则标记
                    CX_plot_event( e, 'dson');
                else
                    CX_plot_event( e, 's'); % 出现标记
                end
            end
        end
    end
    
    %% 统计消失并标记
    if n~=numel(kthdir)
        ind = find(other(:,3)==n);
        if ~isempty(ind)
            for j=ind'
                tmpe = other(j,1);
                if any(divide(:,4)==tmpe) % 如果是母细胞则不算消失
                    continue;
                end
                
                e.x0 = stats(tmpe).Centroid(1);
                e.y0 = stats(tmpe).Centroid(2);
                CX_plot_event( e, 't');
            end
        end
    end
    hold off
    
    %% 保存
    saveas(1,[savepath, strrep(kthdir(n).name,'tif','fig')]);
    close(1);
    
end
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        

            