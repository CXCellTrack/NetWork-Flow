function  CX_plot( e, color )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
alpha1 = e.alpha;
a = e.a;
b = e.b;
x0 = e.x0;
y0 = e.y0;

c=cosd(alpha1);
s=sind(alpha1);
polar_angle=linspace(0,360,361);
xq= a*cosd(polar_angle);
yq= b*sind(polar_angle);
xn=xq*c-yq*s+x0;
yn=xq*s+yq*c+y0;

fill(xn, yn, color,'FaceAlpha', 0.8, 'edgealpha', 0);

end


