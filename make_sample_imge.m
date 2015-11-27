close all

im=imread('C:\Users\Administrator\Desktop\CellTracking论文\merge-split演示图片\t34.png');
imshow(im);
imcontrast(gca); % 调正合适的对比度
set(gca,'position',[0 0 1 1])

% 删除红色椭圆
red_e = findobj('color','r');
delete(red_e);

