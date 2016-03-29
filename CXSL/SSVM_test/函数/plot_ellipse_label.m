function [im, count] = plot_ellipse_label(im, Label, e, scale)

% % 说明是单假说前景，可以用前景形状来代替椭圆 label不能为0，否则变为背景了
ny0 = round(e.y0); ny0(ny0<=0) = 1; ny0(ny0>size(im,1)) = size(im,1);
nx0 = round(e.x0); nx0(nx0<=0) = 1; nx0(nx0>size(im,2)) = size(im,2);
if e.num_hypoth==1 && Label(ny0, nx0)~= 0 
    logi = Label==Label(ny0, nx0);
    count = sum(sum(logi));
    % ------------------------------------------------------------------- %
    % 这个地方的bug在于，有2个十字交叉的区域，椭圆拟合时是独立的，但bwlabel时算作同一类
    % 使得一个椭圆占据2个区域，导致出现
    % The automatic track with label 49 is not consistent with the image data!
    % 因此需要在这添加阈值判断，后面的数据集可以在seglist时消除这种错误
    % 如果区域面积大于椭圆面积太多，则不采用区域表示
    % ------------------------------------------------------------------- %
    if count && count<=1.5*pi*e.a*e.b % count不能为0
        im(logi) = e.color;
        return % 这样做有利于加强分割精度
    end
end

count = 0;
e.a = e.a*scale; % 稍微缩小点椭圆
e.b = e.b*scale;


c=sqrt(e.a^2-e.b^2);   % 焦距
cc=cosd(e.alpha);
ss=sind(e.alpha);
c1=[c*cc+e.x0,c*ss+e.y0];
c2=[-c*cc+e.x0,-c*ss+e.y0];

squares = 5;
for x = floor(e.x0-e.a) : ceil(e.x0+e.a)
    if x<1 || x>size(im,2)
        continue
    end
    for y = floor(e.y0-e.a) : ceil(e.y0+e.a)
        if y<1 || y>size(im,1) % 出界的不算
            continue
        end
        % 要求该点必须不能被别的label占领
        if norm([x,y]-c1) + norm([x,y]-c2) <= 2*e.a && im(y,x)==0
            % 给予椭圆内的点灰度值
            im(y,x) = e.color;
            count = count + 1;
        end
        if count==0 && y==round(e.y0) && x==round(e.x0) % 如果找不到落脚点，则覆盖一个区域，确保能计算AOG和TRA
            im(y-squares:y+squares,x-squares:x+squares) = e.color*ones(squares*2+1);
            count = count + (squares*2+1)^2;
        end
    end
end






