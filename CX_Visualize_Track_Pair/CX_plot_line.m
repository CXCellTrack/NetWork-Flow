function  CX_plot_line( e1, e2, color )

x1 = e1.x0;
y1 = e1.y0;

x2 = e2.x0;
y2 = e2.y0;

if numel(e1.color_index)==1
    thiscolor = color(e1.color_index,:);
else
    thiscolor = [1 1 1];
end
plot([x1 x2], [y1 y2], 'color',thiscolor,'linewidth',1.5);


end