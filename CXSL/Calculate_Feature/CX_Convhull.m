function Area = CX_Convhull(e1, e2)
     
% 先计算e1的圆周点，取361即可
alpha1 = e1.alpha;
a = e1.a;
b = e1.b;
x0 = e1.x0;
y0 = e1.y0;

c=cosd(alpha1);
s=sind(alpha1);
polar_angle=linspace(0,360,50);
xq= a*cosd(polar_angle);
yq= b*sind(polar_angle);
xn1=xq*c-yq*s+x0;
yn1=xq*s+yq*c+y0;

% 再计算e2的圆周点，取20即可
alpha1 = e2.alpha;
a = e2.a;
b = e2.b;
x0 = e2.x0;
y0 = e2.y0;

c=cosd(alpha1);
s=sind(alpha1);
polar_angle=linspace(0,360,50);
xq= a*cosd(polar_angle);
yq= b*sind(polar_angle);
xn2=xq*c-yq*s+x0;
yn2=xq*s+yq*c+y0;

xx = [xn1 xn2];
yy = [yn1 yn2];

[~, Area] = convhull(xx, yy);       
        
end