close all

im=imread('C:\Users\Administrator\Desktop\CellTracking����\merge-split��ʾͼƬ\t34.png');
imshow(im);
imcontrast(gca); % �������ʵĶԱȶ�
set(gca,'position',[0 0 1 1])

% ɾ����ɫ��Բ
red_e = findobj('color','r');
delete(red_e);

