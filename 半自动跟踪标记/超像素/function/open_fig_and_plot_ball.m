function open_fig_and_plot_ball( segpath, fig_dir, bb_dir, t )

[ ~, trackpath ] = getpath( 'competition' );
% ������ͼ�Ĵ�С
screen_size = get(0,'ScreenSize');
first_pos = [1, 50, screen_size(3)/2, screen_size(4)-130];
second_pos = [screen_size(3)/2, 50, screen_size(3)/2, screen_size(4)-130];

% ����ͼƬ
openfig( [ trackpath, '\Pair\���ӻ����ٱ��\', fig_dir(t).name] ); hold on 
if gcf==1
    set(gcf, 'Position', first_pos);
elseif gcf==2
    set(gcf, 'Position', second_pos);
end

set(gca, 'position',[0 0 1 1]);
% ------------------------------------------------------------------- %
% ��ÿ��bsp�ϻ���С��ɫԲ�����ڱ��
bb_path = [ segpath, '\�����ػҶȱ�ǩ\'];
im = imread([bb_path, bb_dir(t).name]); % imshow(im)
% ��ÿ����ǩ�����Ļ�����ɫ���α��
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
