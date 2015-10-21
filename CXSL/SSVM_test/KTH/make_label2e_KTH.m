function stats = make_label2e_KTH( KTH_RES_PATH, e_used, Ellipse )
%
% 将kth中的细胞与我的方法中的椭圆假说一一对应上，以方便计算精度
%

KTH_RES = dir([KTH_RES_PATH, '*.tif']);
stats = cell(numel(KTH_RES),1);

% 如果没有则需要计算一次stats（时间花费较长）
for t=1:numel(KTH_RES)
    im = imread([KTH_RES_PATH, KTH_RES(t).name]);
    erow = e_used(t,e_used(t,:)~=0); % 获取真实用到的椭圆假说

    stats{t} = regionprops(im, 'PixelList');
    for u=1:numel(stats{t})
        stats{t}(u).e = []; % 一个假说点不可能同时位于2个区域，但1个区域内可能有好几个假说点
    end

    disp(['  处理第',num2str(t),'张图片...']);
    tic
    for j=1:numel(erow) % 标记当前区域对应的椭圆假说
        tmpe = erow(j);
        xy = round([Ellipse{t}{tmpe}.x0, Ellipse{t}{tmpe}.y0]);
        % 查找假说中心点是否在区域内
        for u=1:numel(stats{t})
            if isempty(stats{t}(u).PixelList) % 此灰度没出现则跳过
                continue;
            end
            flag = ismember(stats{t}(u).PixelList, xy, 'rows');
            if any(flag==1)
                stats{t}(u).e = [stats{t}(u).e; tmpe];
                break;
            end
        end
    end
    toc
%     num_e = arrayfun(@(x)numel(x.e), stats{t}); % 每个前景中包含的椭圆个数    
end
    
    
    
    
    