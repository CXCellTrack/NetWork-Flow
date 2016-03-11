% ================================================================== %
%
% CX 2015.7.21
% ����ű����ڽ�һ����Ϊ�ӽ�gt��trackdata figure�ֶ�����Ϊground truth figure
% ��Ϊ�ӽ�gt��data������ѵ�����õ���w���ڲ��Լ��ϵĽ��
%
% �������

% ��Ҫ���������1 �����������͵ģ���Ҫɾȥ����켣
%              2 �����޸ĵģ����޸�Ϊ�µ���ȷ�켣
% GUIʹ�÷�����
%
%   1 ѡ����Բ����d����ɾ������Բ����·��ȥ·��ֻ����ǰһ֡�в����������ո��ȷ��
%   2 ѡ�е�һ֡��Բ����ѡ�еڶ�֡��Բ�����һ��Ǩ�ƹ켣�����ո��ȷ��
%   3 ѡ�е�һ֡��Բ�����ո��ȷ���������ʧ�¼�
%   4 ֱ��ѡ�еڶ�֡��Բ�����ո��ȷ������ӳ����¼�
%   5 ѡ�е�һ֡��1������ѡ�еڶ�֡��2�������ո��ȷ�������divide/split�¼�
%   6 ѡ�е�һ֡��2������ѡ�еڶ�֡��1�������ո��ȷ�������merge�¼�
%
%    Attention���������������ʽ�Ĳ��������ᱨ��
% ================================================================== %
clear;close all;
if 0
    dataset = 'competition';
else
    dataset = 'training';
end
[ ~, trackpath ] = getpath( dataset );

fig_addr = [ trackpath, '\Pair\���ӻ����ٱ��\'];
if strcmp(dataset, 'training')
    fig_addr = [ trackpath, '\GT\'];
end
screen_size = get(0,'ScreenSize');
colored_fig_dir = dir([ fig_addr, '*.fig' ]);  
frame = numel(colored_fig_dir);


% ����ȫ�ֱ��� GT_move �洢��׼Ǩ�ƴ�
% ����ȫ�ֱ��� GT_delete �洢��ǰ�������ҪҪ��ɾȥ���¼�
global GT_move GT_delete global_x t;

mkdir([ trackpath, '\GT']);
handGt = [ trackpath, '\GT\Hand_GT_New.mat'];
if exist(handGt, 'file')
    load( handGt );
    GT_move = GT_move_s;
    GT_delete = GT_delete_s;
else
    GT_move = cell(frame-1,1);
    GT_delete = cell(frame-1,1);
    
end


for t = 1:frame-1 % 2015.10.4��hela1ѵ������ǰ80֡������

    % �������2������ֹ�����޸�ʱ����
    GT_move{t} = zeros(50,4); % Ԥ��20��Ӧ�ù��õ�
    GT_delete{t} = [];
    
    % ��2������ͼƬ
    for ss = t:t+1
        fig_name = [ fig_addr, colored_fig_dir(ss).name ];
        openfig(fig_name, 'new', 'visible');
        set(gcf, 'Position', screen_size);
    end
    % ȡ��figure����Բ�ľ��e1,e2
    h1 = get(1, 'children');
    h2 = get(2, 'children');
    e1 = findobj(h1, 'Type', 'Line');
    e2 = findobj(h2, 'Type', 'Line');

    % ���ð�����Ӧ�͵����Ӧ
    global_x = 1;
    set(1, 'keypressfcn', @pressBlank); % ���ո��������
    set(2, 'keypressfcn', @pressBlank);
    set(e1, 'ButtonDownFcn', @myclick_1);
    set(e2, 'ButtonDownFcn', @myclick_2);

    disp(['  ���ڱ��',num2str(t),'��',num2str(t+1),'֡']);
    
    % ͨ��input���������������ʱ��
    keyin = input('  �ڿ���̨�ϰ�enter����������֡\n', 's');
    if strcmp(keyin, '')
        close(1);
        close(2);
    end
    
    if 1
        disp('  ��������Hand_GT_New.mat ��');
        GT_delete_s = GT_delete;
        GT_move_s = GT_move;
        save([ trackpath, '\GT\Hand_GT_New.mat'], 'GT_move_s','GT_delete_s');
    end

end






















