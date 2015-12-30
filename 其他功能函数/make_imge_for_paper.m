close all

im=imread('C:\Users\Administrator\Desktop\CellTracking论文\merge-split演示图片\t34.png');
imshow(im);
imcontrast(gca); % 调正合适的对比度
set(gca,'position',[0 0 1 1])

% 删除红色椭圆
red_e = findobj('color','r');
delete(red_e);

%% 在hela1中找到PSL和OURS-P错误的地方

im=imread('C:\Users\Administrator\Desktop\PSL错误\OURS\t80.tif');
% im = mat2gray(im);
% h = fspecial('gaussian', 10); 
% im = imfilter(im, h);
imshow(im)
% imshow(imadjust(im));
imshow(mat2gray(im));


im = zeros(700,1100);imshow(im);hold on
load('E:\datasets\first_edition\competition_datasets\N2DL-HeLa\01_2-16_track\Pair\Pre_data_New.mat');
remain = [224,225,226,210,240,245,237];
t = 81;
rng(1);
for j=remain
    CX_fill_color( Ellipse{t}{j}, hsv, randi(64) );
end
hold off












