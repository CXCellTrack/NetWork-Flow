clear

% 注意到第一届的6个sim中，seg的答案都是全的，因此可直接使用标准答案进行
simpath = 'E:\datasets\first_edition\training_datasets\N2DH-SIM\';
for ii=1:6
    gtpath = sprintf('%s0%d_GT\\TRA\\',simpath,ii); % sprintf遇到\需要转义
    segpath = sprintf('%s0%d_0-00_seg\\',simpath,ii);
    gtdir = dir([gtpath,'*.tif']);
    for jj=1:numel(gtdir)
        im = imread([gtpath, gtdir(jj).name]);
        im(im~=0) = 1; im = logical(im); % 将gt中的灰度图转换为二值图像
        imwrite(im, [segpath, gtdir(jj).name]);
    end
end

