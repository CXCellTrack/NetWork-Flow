function eval_str = CX_plot_event_sp( eval_str, sp ,char)
%
% ����im�����䣬��ʼshow�Ļ���Щ��ǻ���ȥ�޷��ƶ���������ͼ�� 
% ��˲������ְ취���ȼ�¼�»�ͼ��䣬��im�������show��Ȼ��һ�������
% 2015.12.14

x0 = sp.centroid(1);
y0 = sp.centroid(2);

switch char  
    case 'd' % ĸϸ������Ҫdivide
        eval_str = [ eval_str, 'plot(',num2str(x0),',',num2str(y0),',''ws'',''LineWidth'', 1);' ];
    case 'dson' % ��ϸ��
        eval_str = [ eval_str, 'plot(',num2str(x0),',',num2str(y0),',''ks'',''LineWidth'', 1);' ];
    case 's' % �³��ֵ�
        eval_str = [ eval_str, 'plot(',num2str(x0),',',num2str(y0),',''wo'',''LineWidth'', 1);' ];
    case 't' % ��һ֡��ʧ
        eval_str = [ eval_str, 'plot(',num2str(x0),',',num2str(y0),',''wx'',''LineWidth'', 1);' ];
    case 'm' % �ں϶���
        eval_str = [ eval_str, 'plot(',num2str(x0),',',num2str(y0),',''r+'',''LineWidth'', 1);' ];
    case 'mson' % �ں϶���
        eval_str = [ eval_str, 'plot(',num2str(x0),',',num2str(y0),',''k+'',''LineWidth'', 1);' ]; 
    case 'v' % ��Ҫsplit
        eval_str = [ eval_str, 'plot(',num2str(x0),',',num2str(y0),',''rd'',''LineWidth'', 1);' ]; 
    case 'vson' % split������Сϸ��
        eval_str = [ eval_str, 'plot(',num2str(x0),',',num2str(y0),',''kd'',''LineWidth'', 1);' ]; 
        
end

end

