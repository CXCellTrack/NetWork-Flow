close all

im=imread('C:\Users\Administrator\Desktop\CellTracking����\merge-split��ʾͼƬ\t34.png');
imshow(im);
imcontrast(gca); % �������ʵĶԱȶ�
set(gca,'position',[0 0 1 1])

% ɾ����ɫ��Բ
red_e = findobj('color','r');
delete(red_e);

%% ��hela1���ҵ�PSL��OURS-P����ĵط�

im=imread('C:\Users\Administrator\Desktop\PSL����\OURS\t80.tif');
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












