function CX_plot_event( e ,char)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
x0 = e.x0;
y0 = e.y0;
switch char  
    case 'd' % ĸϸ������Ҫdivide
        plot(x0, y0,'ws','LineWidth', 1);
    case 'dson' % ��ϸ��
        plot(x0, y0,'ks','LineWidth', 1);
    case 's' % �³��ֵ�
        plot(x0, y0,'wo','LineWidth', 1);
    case 't' % ��һ֡��ʧ
        plot(x0, y0,'wx','LineWidth', 1);
    case 'm' % �ں϶���
        plot(x0, y0,'r+','LineWidth', 1);
    case 'mson' % �ں϶���
        plot(x0, y0,'k+','LineWidth', 1);   
    case 'v' % ��Ҫsplit
        plot(x0, y0,'rd','LineWidth', 1);
    case 'vson' % split������Сϸ��
        plot(x0, y0,'kd','LineWidth', 1);
        
end

end

