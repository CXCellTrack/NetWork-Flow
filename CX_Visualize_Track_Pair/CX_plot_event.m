function CX_plot_event( e ,char)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
x0 = e.x0;
y0 = e.y0;
switch char  
    case 'd'
        plot(x0, y0,'ws','LineWidth', 1);
    case 's'
        plot(x0, y0,'wo','LineWidth', 1);
    case 't'
        plot(x0, y0,'wx','LineWidth', 1);
    case 'm'
        plot(x0, y0,'b+','LineWidth', 1);
    case 'v'
        plot(x0, y0,'bd','LineWidth', 1);
end

end

