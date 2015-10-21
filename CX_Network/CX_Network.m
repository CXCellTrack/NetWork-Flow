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
tolerance = 2.0;
remove_small = 25; % �Ƴ�С�������ֵ������΢��һ�㣡

for frame=1:length(rawpic_dir)
    
    iteration_num = -1;
    lunkuo = 'error';
    pic = [ segpath, '\', rawpic_dir(frame).name ];

    while strcmp(lunkuo, 'error') % ��������������˲�����
        iteration_num = iteration_num + 1;
        if iteration_num>=4
            error('�˲������ﵽ4�Σ��鿴�������');
        end
%         iteration_num = 1;
        [ ellipse{frame}, lunkuo] = CX_Fit( frame, pic, iteration_num, tolerance, remove_small);   %%�������˲�
    end

    % ����FOI�ڵ�������������Ҫreplot
    lunkuo_name = [ lunkuo_addr, rawpic_dir(frame).name ];
    imwrite(lunkuo, lunkuo_name);
    
    % ������Ϻ��figure
    disp('������');
    savename = strcat(output_addr, rawpic_dir(frame).name(1:end-4), '_fit.fig');
    saveas(1, savename);  % �������
    close('1');
    save(ellipse_address, 'ellipse');
end
