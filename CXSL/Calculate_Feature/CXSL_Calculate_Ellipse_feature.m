function  e = CXSL_Calculate_Ellipse_feature( e, dataset_addr )
% structured learing部分——特征提取
%   Detailed explanation goes here
frame = numel(e);
raw_dir = dir([ dataset_addr, '\*.tif' ]);

for t=1:frame
    raw_im = imread([ dataset_addr, '\', num2str(raw_dir(t).name) ]);
    raw_im_db = mat2gray(raw_im);
    [height,~]=size(e{t});

    for n=1:height        %%每个椭圆 此处处理的是拉直过后的cell
        % 问题椭圆在数据预处理中已经滤去
%             if e{t}{n}.status == 0  %% 椭圆有问题则跳过
%                 continue;
%             end
%             e_4_c = struct2array(e{t}{n}); % 转化为array，供c语言调用
%             目前有一点问题，无法循环完所有的t
%             e{t}{n}.feature.intensity = CXSL_feature_C(e_4_c, raw_im_db);
        e{t}{n}.feature.intensity = Cal_Intensity_feature(e{t}{n}, raw_im_db);

        % 椭圆基本信息
        e{t}{n}.feature.geo.position = [ e{t}{n}.x0, e{t}{n}.y0 ]; % 位置
        e{t}{n}.feature.geo.size = pi * e{t}{n}.a * e{t}{n}.b;     % 尺寸
        e{t}{n}.feature.geo.eccentric = sqrt( e{t}{n}.a^2 - e{t}{n}.b^2 )/e{t}{n}.a;  % 偏心率
        e{t}{n}.feature.geo.alpha = e{t}{n}.alpha;                   % 偏角
        e{t}{n}.feature.dist2border = min([size(raw_im_db)-[e{t}{n}.y0, e{t}{n}.x0], e{t}{n}.y0, e{t}{n}.x0]);
        
    end % end for
end % end frame

end
    
            

%######################### matlab中 1<=a<=2 表示的不是 1<=a && a<=2  注意！！！

function [ intensity ] = Cal_Intensity_feature(e1, im)
%
%
bins = 16; % 直方图的通道数目
[h, w] = size(im);
% pixellist 用于记录椭圆内的点灰度值（以后也可以连坐标一起记录！2015.6.29）
pixellist = [];

c1=sqrt(e1.a^2-e1.b^2);   %%焦距
cc1=cosd(e1.alpha);
ss1=sind(e1.alpha);
c11=[c1*cc1+e1.x0,c1*ss1+e1.y0];
c12=[-c1*cc1+e1.x0,-c1*ss1+e1.y0];
count=0;
for x=e1.x0-e1.a : e1.x0+e1.a
    for y=e1.y0-e1.a : e1.y0+e1.a
        xx = round(x);
        yy = round(y);
        if yy>h || xx>w || yy<1 || xx<1 % 在边界的椭圆，找点可能会找到图外面去
            continue;
        end
        
        if sqrt((x-c11(1))^2+(y-c11(2))^2)+sqrt((x-c12(1))^2+(y-c12(2))^2)<=2*e1.a
            % 保存椭圆内的点灰度值
            count = count + 1;
            pixellist(count) = im(yy,xx);
        end
    end
end

intensity.hist = hist(pixellist, bins)'/count;
intensity.sum = sum(pixellist);
intensity.mean = mean(pixellist);
intensity.devia = std(pixellist);

end

            




