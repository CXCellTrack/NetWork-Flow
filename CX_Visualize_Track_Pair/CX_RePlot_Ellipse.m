function  CX_RePlot_Ellipse( dataset )

% ----------------------------------------------------------------------- %
% 这个函数用于重新绘制优化后椭圆，去掉了分段直线和优化时去掉的椭圆
% 最终在“新拟合图”这个文件夹中生成新椭圆图
% 绘制可视化跟踪的图像就是在新椭圆图的基础上
% ----------------------------------------------------------------------- %

% if 1
%     dataset = 'competition'; % 选择训练还是测试
% else
%     dataset = 'training';
% end
[ segpath trackpath ] = getpath( dataset );
% 载入椭圆数据

load([trackpath, '\Pair\Pre_data_new.mat'], 'Ellipse');
lunkuo_addr = [segpath, '\FOI提取轮廓\'];  % 只需要修改此处
lunkuo_dir = dir([ lunkuo_addr, '*.tif' ]);
output_addr = [trackpath, '\新拟合图\'];

frame = numel(Ellipse);

for t=1:frame
    disp(['处理第',num2str(t),'帧...']);
    pic_name = [ lunkuo_addr, lunkuo_dir(t).name ];
    edgeim = imread(pic_name);
    % 读入轮廓图片作为底板
    figure;
    imshow(edgeim);hold;
    % 绘制椭圆
    for j=1:numel(Ellipse{t})
        e = Ellipse{t}{j};
        alpha1 = e.alpha;
        a = e.a;
        b = e.b;
        x0 = e.x0;
        y0 = e.y0;

        c = cosd(alpha1);
        s = sind(alpha1);
        polar_angle = linspace(0,360,361);
        xq = a*cosd(polar_angle);
        yq = b*sind(polar_angle);
        xn = xq*c-yq*s+x0;
        yn = xq*s+yq*c+y0;
        plot(xn,yn,'g','LineWidth',1.5,'displayname',num2str(j));
        
    end
    hold;
    savename= strcat( output_addr, lunkuo_dir(t).name(1:end-4), '_fit.fig');
    saveas(1,savename);
    close('1');
end

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    