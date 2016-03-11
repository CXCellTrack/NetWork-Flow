clear

% ע�⵽��һ���6��sim�У�seg�Ĵ𰸶���ȫ�ģ���˿�ֱ��ʹ�ñ�׼�𰸽���
simpath = 'E:\datasets\first_edition\training_datasets\N2DH-SIM\';
for ii=1:6
    gtpath = sprintf('%s0%d_GT\\TRA\\',simpath,ii); % sprintf����\��Ҫת��
    segpath = sprintf('%s0%d_0-00_seg\\',simpath,ii);
    gtdir = dir([gtpath,'*.tif']);
    for jj=1:numel(gtdir)
        im = imread([gtpath, gtdir(jj).name]);
        im(im~=0) = 1; im = logical(im); % ��gt�еĻҶ�ͼת��Ϊ��ֵͼ��
        imwrite(im, [segpath, gtdir(jj).name]);
    end
end

