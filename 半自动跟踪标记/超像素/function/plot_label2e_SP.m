function plot_label2e_SP( stats, dataset )


[ segpath, trackpath ] = getpath( dataset );
bb_path = [ segpath, '\�����ػҶȱ�ǩ\'];
bb_dir = dir([bb_path, '\*.png']);

for t=1:numel(bb_dir)
    disp(['  ����GT��ǩͼ', num2str(t), '...']);
    im = imread([bb_path, bb_dir(t).name]); % imshow(im)
    bb = zeros(size(im)); % ȫ�ڻ�ͼ�װ�
    col = 1;
    maskim = drawregionboundaries(im, bb, col);
    imshow(maskim); % ���Ʊ�ǩͼ
    hold on
    
    % ��ÿ����ǩ�����Ļ�����ɫ���α��
    stats_SP = regionprops(im, 'centroid');
    for j=1:numel(stats_SP)
        x0 = stats_SP(j).Centroid(1);
        y0 = stats_SP(j).Centroid(2);
        rr = 10;
        polar_angle = linspace(0,360,361);
        xq = rr*cosd(polar_angle) + x0;
        yq = rr*sind(polar_angle) + y0;
        plot(xq, yq, 'g','DisplayName',num2str(j),'linewidth',1.5);
    end
    
    if strcmp(dataset,'training');
        % ��gtÿ����ǩ�����Ļ��ϰ�ɫ*��
        for label=1:numel(stats{t})
            if isnan( stats{t}(label).Centroid(1) ) % �������ΪNaN��������
                continue;   
            end
            x = stats{t}(label).Centroid(1);
            y = stats{t}(label).Centroid(2);
            plot(x, y, 'w*', 'DisplayName',num2str(label));
    %         text(x, y, num2str(label), 'color', [1 1 1], 'fontsize', 8);
        end
    end
    
    hold off 
    
    % ����ͼƬ
    saveas(1, [ trackpath, '\GT\label_and_e\', strrep(bb_dir(t).name,'png','fig')]);
    close('1');
    
end














