function CX_plot_event( e ,char)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
x0 = e.x0;
y0 = e.y0;
switch char  
    case 'd' % 母细胞，将要divide
        plot(x0, y0,'ws','LineWidth', 1);
    case 'dson' % 子细胞
        plot(x0, y0,'ks','LineWidth', 1);
    case 's' % 新出现的
        plot(x0, y0,'wo','LineWidth', 1);
    case 't' % 下一帧消失
        plot(x0, y0,'wx','LineWidth', 1);
    case 'm' % 融合而来
        plot(x0, y0,'r+','LineWidth', 1);
    case 'mson' % 融合而来
        plot(x0, y0,'k+','LineWidth', 1);   
    case 'v' % 将要split
        plot(x0, y0,'rd','LineWidth', 1);
    case 'vson' % split而来的小细胞
        plot(x0, y0,'kd','LineWidth', 1);
        
end

end

