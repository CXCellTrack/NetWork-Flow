


I = imadjust(im2double(imread('E:\datasets\first editon\competition datasets\N2DL-HeLa\01\t70.tif')));
w=466; h=382;
w1=228; h1=632;


[J,reg_size,reg_mean]=regiongrowing(I,h,w,2^13/(2^16-1));
[J1,reg_size1,reg_mean1]=regiongrowing(I,h1,w1,2^10/(2^16-1));imshow(I+J+J1)

kmeans


figure, imshow(I+J);

