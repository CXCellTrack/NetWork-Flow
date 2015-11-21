function  CX_plot_ball( e, color )

x0 = e.x0;
y0 = e.y0;

rr = 8;
polar_angle=linspace(0,360,361);
xq= rr*cosd(polar_angle) + x0;
yq= rr*sind(polar_angle) + y0;

if numel(color)==3
    thiscolor = color;
else
    thiscolor = color(e.color_index,:);
end
fill(xq, yq, thiscolor,'FaceAlpha', 1, 'edgealpha', 0);

end