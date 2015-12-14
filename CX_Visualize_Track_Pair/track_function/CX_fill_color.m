function  CX_fill_color( e, color, index )

% 往椭圆中填充颜色
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

if ~exist('index', 'var')
    thiscolor = color; % 也可以直接输入color不输入index
else  
    assert(numel(index)==1)
    thiscolor = color(index,:);
end

% plot(xn, yn, 'g','linewidth',1.5);
fill(xn, yn, thiscolor,'FaceAlpha', 1, 'edgealpha', 0); %, 'linewidth',1.5, 'edgecolor','g'

end


