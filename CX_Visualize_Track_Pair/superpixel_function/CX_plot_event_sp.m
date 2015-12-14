function eval_str = CX_plot_event_sp( eval_str, sp ,char)
%
% 由于im后面会变，开始show的话这些标记画上去无法移动到后来的图上 
% 因此采用这种办法，先记录下画图语句，等im不变后再show，然后一起做标记
% 2015.12.14

x0 = sp.centroid(1);
y0 = sp.centroid(2);

switch char  
    case 'd' % 母细胞，将要divide
        eval_str = [ eval_str, 'plot(',num2str(x0),',',num2str(y0),',''ws'',''LineWidth'', 1);' ];
    case 'dson' % 子细胞
        eval_str = [ eval_str, 'plot(',num2str(x0),',',num2str(y0),',''ks'',''LineWidth'', 1);' ];
    case 's' % 新出现的
        eval_str = [ eval_str, 'plot(',num2str(x0),',',num2str(y0),',''wo'',''LineWidth'', 1);' ];
    case 't' % 下一帧消失
        eval_str = [ eval_str, 'plot(',num2str(x0),',',num2str(y0),',''wx'',''LineWidth'', 1);' ];
    case 'm' % 融合而来
        eval_str = [ eval_str, 'plot(',num2str(x0),',',num2str(y0),',''r+'',''LineWidth'', 1);' ];
    case 'mson' % 融合而来
        eval_str = [ eval_str, 'plot(',num2str(x0),',',num2str(y0),',''k+'',''LineWidth'', 1);' ]; 
    case 'v' % 将要split
        eval_str = [ eval_str, 'plot(',num2str(x0),',',num2str(y0),',''rd'',''LineWidth'', 1);' ]; 
    case 'vson' % split而来的小细胞
        eval_str = [ eval_str, 'plot(',num2str(x0),',',num2str(y0),',''kd'',''LineWidth'', 1);' ]; 
        
end

end

