function [im, count] = plot_ellipse_label(im, Label, e, scale)

% % 说明是单假说前景，可以用前景形状来代替椭圆 label不能为0，否则变为背景了
if e.num_hypoth==1 && Label(round(e.y0),round(e.x0))~= 0 
    logi = Label==Label(round(e.y0),round(e.x0));
    count = sum(sum(logi));
    % ------------------------------------------------------------------- %
    % 这个地方的bug在于，有2个十字交叉的区域，椭圆拟合时是独立的，但bwlabel时算作同一类
    % 使得一个椭圆占据2个区域，导致出现
    % The automatic track with label 49 is not consistent with the image data!
    % 因此需要在这添加阈值判断，后面的数据集可以在seglist时消除这种错误
    % 如果区域面积大于椭圆面积太多，则不采用区域表示
    % ------------------------------------------------------------------- %
    if count<=1.5*pi*e.a*e.b 
        im(logi) = e.color;
        return % 这样做有利于加强分割精度
    end
end

count = 0;
e.a = e.a*scale; % 稍微缩小点椭圆
e.b = e.b*scale;

for x = floor(e.x0-e.a) : ceil(e.x0+e.a)
    if x<1 || x>size(im,2)
        continue
    end
    for y = floor(e.y0-e.a) : ceil(e.y0+e.a)
        if y<1 || y>size(im,1) % 出界的不算
            continue
        end
        
        c=sqrt(e.a^2-e.b^2);   % 焦距
        cc=cosd(e.alpha);
        ss=sind(e.alpha);
        c1=[c*cc+e.x0,c*ss+e.y0];
        c2=[-c*cc+e.x0,-c*ss+e.y0];
        
        if norm([x,y]-c1) + norm([x,y]-c2) <= 2*e.a
            % 给予椭圆内的点灰度值
            im(y,x) = e.color;
            count = count + 1;
        end
    end
end

