% ͳ�Ʒָ��

clear
% �����׼�𰸵����ĵ�����
load('E:\datasets\first_edition\training_datasets\N2DH-SIM\05_0-00_track\GT\center_gt.mat');
load('E:\datasets\first_edition\training_datasets\N2DH-SIM\05_0-00_track\Pair\Pre_data_new.mat','n');
contour_path = 'E:\datasets\first_edition\training_datasets\N2DH-SIM\05_0-00_seg\FOI��ȡ����\';

frame = numel(n);
contour_dir = dir([ contour_path, '*.tif' ]);

%% ����һһ��Ӧ����Ӧ���ϵĿɼ���prec��recall
center_e = cell(frame,1);
distance = cell(frame,1);
TP = zeros(frame,1);   % ǰ����ϸ����Ӧ����
Tnum = zeros(frame,1); % ʵ��ϸ������
Pnum = zeros(frame,1); % �ָ�õ���ǰ������

for t=1:frame % ��frame���У���Ŀǰ�����û�еĽ��б��
    disp(['�����',num2str(t),'֡...']);
    center_gt{t} = center_gt{t}(~isnan(center_gt{t}(:,1)),:); % ȥ��center��Ϊnan����Щ��
    
    im = imread([contour_path, contour_dir(t).name]); % imshow(im)
    im = bwareaopen(im, 50);
    stats = regionprops(im, 'centroid','MajorAxisLength');
    center_region = cat(1, stats.Centroid);
    % Ϊÿ��ǰ���ҵ����������*�㣬��ôÿ��ǰ��������һ��*�㣬��ʱΪTP�����Ҳ���*�㣬ΪFP���������*��ΪFN
    distance{t} = dist(center_region, center_gt{t}'); 
    for re=1:numel(stats)
        [thisD, label] = min(distance{t}(re,:)); 
        if thisD<=stats(re).MajorAxisLength % �������С������Բ�ģ���*��λ�������ڵĿ��ٽ����жϣ�
            TP(t) = TP(t) + 1;
        end
    end
    
    [ Pnum(t), Tnum(t) ] = size(distance{t});
    
end
     
precision = sum(TP)/sum(Pnum);
recall = sum(TP)/sum(Tnum);
F_m = precision*recall*2/(precision+recall);
% ��ӡ���
fprintf('precision:\t%f\nrecall:\t\t%f\nF_m:\t\t%f\n',precision,recall,F_m);
    
    
    
    
    
    
    
    
    
    
    
    
    
    