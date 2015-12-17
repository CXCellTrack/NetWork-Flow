function open_fig_and_plot_ball( segpath, fig_dir, bb_dir, t )

[ ~, trackpath ] = getpath( 'competition' );
% 调正幅图的大小
screen_size = get(0,'ScreenSize');
first_pos = [1, 50, screen_size(3)/2, screen_size(4)-130];
second_pos = [screen_size(3)/2, 50, screen_size(3)/2, screen_size(4)-130];

% 载入图片
openfig( [ trackpath, '\Pair\可视化跟踪标记\', fig_dir(t).name] ); hold on 
if gcf==1
    set(gcf, 'Position', first_pos);
elseif gcf==2
    set(gcf, 'Position', second_pos);
end

set(gca, 'position',[0 0 1 1]);
% ------------------------------------------------------------------- %
% 在每个bsp上画上小绿色圆，用于标记
bb_path = [ segpath, '\超像素灰度标签\'];
im = imread([bb_path, bb_dir(t).name]); % imshow(im)
% 在每个标签的中心画上绿色方形标记
stats_SP = regionprops(im, 'centroid');
for j=1:numel(stats_SP)
    x0 = stats_SP(j).Centroid(1);
    y0 = stats_SP(j).Centroid(2);
    rr = 10;
    polar_angle = linspace(0,360,361);
    xq = rr*cosd(polar_angle) + x0;
    yq = rr*sind(polar_angle) + y0;
    plot(xq, yq, 'k','DisplayName',num2str(j),'linewidth',1.5);
end

hold off
% ------------------------------------------------------------------- %
