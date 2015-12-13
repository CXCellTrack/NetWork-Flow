function plot_label2e( stats, fig_dir )

[ ~, trackpath ] = getpath( 'training' );
% ��Ŀ���ļ����ڵ�ͼƬ��������frameʱ�����л��Ʋ�����
for t=1:numel(stats)
    disp(['����GTͼƬ', num2str(t), '...']);
    openfig([ trackpath, '\�����ͼ\', fig_dir(t).name ]);
    hold;   

    for label=1:numel(stats{t})
        if isnan( stats{t}(label).Centroid(1) ) % �������ΪNaN��������
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
% �򿪸����ɵ�GTͼƬ���鿴��Ե���Ƿ��ܹ���Ӧ��
if 0
    output_dir = dir([ trackpath, '\GT\*.fig']);
    for t=51:frame
        openfig([ trackpath, '\GT\', output_dir(t).name]);
    end
end