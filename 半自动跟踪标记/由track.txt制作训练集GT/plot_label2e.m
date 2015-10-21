function plot_label2e( stats, fig_dir )

[ ~, trackpath ] = getpath( 'training' );
% 当目标文件夹内的图片数量不到frame时，进行绘制并保存
for t=1:numel(stats)
    disp(['绘制GT图片', num2str(t), '...']);
    openfig([ trackpath, '\新拟合图\', fig_dir(t).name ]);
    hold;   

    for label=1:numel(stats{t})
        if isnan( stats{t}(label).Centroid(1) ) % 如果该列为NaN，则跳过
            continue;   
        end
        x = stats{t}(label).Centroid(1);
        y = stats{t}(label).Centroid(2);
        plot(x, y, 'w*', 'DisplayName',num2str(label));
%         text(x, y, num2str(label), 'color', [1 1 1], 'fontsize', 8);
    end
    hold;
    saveas(1, [ trackpath, '\GT\label_and_e\', fig_dir(t).name]);
    close('1');
end
    % #################################################

    
    
% ---------------------------------------- %
% 打开刚生成的GT图片，查看边缘处是否能够对应上
if 0
    output_dir = dir([ trackpath, '\GT\*.fig']);
    for t=51:frame
        openfig([ trackpath, '\GT\', output_dir(t).name]);
    end
end