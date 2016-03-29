%######################################
%
% 2015.5.30 CX on desk
% ���ã�����ű����ڴӳ�ʼ�ָ��еõ���Բ��ϼ�˵
% ���ݴ洢����Բ���ݱ����� raw_ellipse ��
% ������ϵ�� �������� CX_fit ����ÿ֡ͼƬ����Ͻ�� 
% 2015.6.5�������� FOI ��ȡ��ʹ��Ե���жϷ��ϱ�׼ ��
%
%######################################
clear;close all;

if 1
    dataset = 'competition'; % ѡ��ѵ�����ǲ���
else
    dataset = 'training';
end
[ segpath, ~ ] = getpath( dataset );

segpath = [segpath, '_gt'];

rawpic_dir=dir([ segpath, '\*.tif' ]); % ԭʼtifͼƬ��ַ
% �����Ϻ�ͼƬ��ַ 
output_addr = [ segpath, '\FOI���ͼ2.0\'];
if ~exist(output_addr, 'dir')
    mkdir(output_addr);
end
    
% ���FOI��ͼƬ�����ĵ�ַ
lunkuo_addr=[ segpath, '\FOI��ȡ����\'];
if ~exist(lunkuo_addr, 'dir')
    mkdir(lunkuo_addr);
end
    

% ���ⲿ����e
ellipse_address = [ output_addr, 'raw_ellipse.mat' ];
if exist(ellipse_address, 'file')
    load(ellipse_address);
else
    ellipse = cell(length(rawpic_dir),1);
end
%########## ��Ϲ��� �Ƴ�С���� 
tolerance = 5.0; % �򵥵����ݼ��ɽ�tol����
remove_small = 20; % �Ƴ�С�������ֵ������΢��һ�㣡

%% ��ʼѭ������ͼƬ 

for t=24:length(rawpic_dir)
    
    iteration_num = -1;
    lunkuo = 'error';
    pic = [ segpath, '\', rawpic_dir(t).name ];

    while strcmp(lunkuo, 'error') % ��������������˲�����
        iteration_num = iteration_num + 1;
        if iteration_num>5
            error('�˲������ﵽ5�Σ��鿴�������');
        end
%         iteration_num = 1;
        [ ellipse{t}, lunkuo] = CX_Fit( t, pic, iteration_num, tolerance, remove_small); % �������˲�
    end

    % ����FOI�ڵ�������������Ҫreplot
    lunkuo_name = [ lunkuo_addr, rawpic_dir(t).name ];
    imwrite(lunkuo, lunkuo_name);
    
    % ������Ϻ��figure
    disp('������');
    savename = strcat(output_addr, rawpic_dir(t).name(1:end-4), '_fit.fig');
    saveas(1, savename);  % �������
    close('1');
    save(ellipse_address, 'ellipse');
end

% end
