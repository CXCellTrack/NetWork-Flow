function  SuperPixel = CXSL_Calculate_SP_feature( SuperPixel, segpath, n )
% structured learing部分——特征提取
% 提取单个Superpixel所包含的特征 
%

% 需要给出原始图片的地址
last = max(strfind(segpath, '\'));
raw_addr = segpath(1:last+2);  
raw_dir = dir([ raw_addr, '\*.tif' ]);

frame = numel(SuperPixel);

%% 开始循环：计算每个SP的特征
for t=1:frame % 每帧
    if isempty(SuperPixel{t}) % 滤去空cell
        continue;
    end
    disp(['  计算第',num2str(t),'帧...']);
    im = imread([ raw_addr, '\', num2str(raw_dir(t).name) ]);
    im_db = mat2gray(im);

    for j=1:n(t) % 每个SP 
        SP = SuperPixel{t}{j};
        % baseboard 绘图底板，在底板上画出SP的位置
        bb = zeros(size(im));
        for pl=1:size(SP.pixellist,1)
            bb(SP.pixellist(pl,2), SP.pixellist(pl,1)) = 1;
        end % imshow(bb)
        % 计算SP的特征
        spfeature = regionprops(bb, im_db, 'MeanIntensity','PixelValues',... % 灰度特征
            'Area','Centroid',... % 面积、质心
            'Eccentricity','MajorAxisLength','MinorAxisLength','Orientation',... % 类椭圆特征
            'Extent','Solidity');
        spfeature.DeviaIntensity = std(spfeature.PixelValues); % 再添加灰度标准差
        spfeature.SumIntensity = sum(spfeature.PixelValues); % 灰度和
        spfeature.HistIntensity = hist(spfeature.PixelValues, 16)'/spfeature.Area; % 灰度直方图
        % 到边距离
        spfeature.Dist2border = min([ size(im_db)-[SP.centroid(2), SP.centroid(1)], SP.centroid(2), SP.centroid(1) ]);
        % 移除像素值
        spfeature = rmfield(spfeature,'PixelValues');
        % 赋值，删除其他属性，只保留特征
        SuperPixel{t}{j} = spfeature;
        
    end
end % end frame

    
            



