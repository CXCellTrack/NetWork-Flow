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

if 0
    dataset = 'competition'; % ѡ��ѵ�����ǲ���
else
    dataset = 'training';
end
[ segpath, ~ ] = getpath( dataset );

raw_dir = dir([segpath(1:end-9), '\*.tif']);
bw_dir = dir([ segpath, '\*.tif' ]); % ��ʼ�ָ�tifͼƬ��ַ

% ���������ͼƬ��ַ 
output_addr1 = [ segpath, '\�����ػҶȱ�ǩ\'];
if ~exist(output_addr1, 'dir')
    mkdir(output_addr1);
end
output_addr2 = [ segpath, '\�����ز�ɫ��ǩ\'];
if ~exist(output_addr2, 'dir')
    mkdir(output_addr2);
end
output_addr3 = [ segpath, '\����������ͼ\'];
if ~exist(output_addr3, 'dir')
    mkdir(output_addr3);
end

% ���ⲿ����SuperPixel.mat
SP_addr = [ output_addr1, 'SuperPixel.mat' ];
if exist(SP_addr, 'file')
    load(SP_addr);
else
    SuperPixel = cell(length(raw_dir),1);
end
%% ��ʼѭ������ͼƬ 
nsp = 200; % 200��������

for frame=1:6 % length(bw_dir)
    tic
    % ����ͼƬ��ַ
    raw_pic = [segpath(1:end-9), '\', raw_dir(frame).name];
    bw_pic = [segpath, '\', bw_dir(frame).name];
    disp(['  ����',raw_pic ]);
    % ���ɱ�ǩͼƬ�ͳ�����cell
    [ SP, new_labels, maskim, RGB_label ] = CX_Generate_SP_in_1( raw_pic, bw_pic, nsp );
    
    % ����õ���ͼƬ
    disp('  ������');
    savename1 = [ output_addr1, bw_dir(frame).name(1:end-4), '_sp.png' ];
    imwrite(uint16(new_labels), savename1);
    savename2 = [ output_addr2, bw_dir(frame).name(1:end-4), '_sp_color.png' ];
    imwrite(RGB_label, savename2);
    savename3 = [ output_addr3, bw_dir(frame).name(1:end-4), '_sp.png' ];
    imwrite(maskim, savename3);
    % ���泬���ؼ�˵����
    SuperPixel{frame} = SP;
    save(SP_addr, 'SuperPixel');
    toc
end








