% ================================================================== %
%
% CX 2015.12.13
% ����ű����ڴ�ԭʼ�ĳ����ؼ�˵�б�ǳ���ȷ��
%
% �������

% ��Ҫ���������1 �����������͵ģ���Ҫɾȥ����켣
%              2 �����޸ĵģ����޸�Ϊ�µ���ȷ�켣
% GUIʹ�÷�����
%
%   1 Ǩ�ƣ���һ֡�е������bsp���csp���ڶ�֡�е������bsp���csp�����ո��ȷ��
%   2 ��ʧ����һ֡�е������bsp���csp�����ո��ȷ��
%   3 ���֣��ڶ�֡�е������bsp���csp�����ո��ȷ��
%   4 ���ѣ���һ֡�е������bsp���cspĸϸ�����ڶ�֡ѡ���һ����ϸ����Ҫ��cȷ������ѡ�ڶ�������󰴿ո��ȷ��
%   5 �ϲ�����һ֡ѡ���һ��Դϸ����Ҫ��cȷ������ѡ�ڶ������ڶ�֡�е������bsp���csp��ϸ������󰴿ո��ȷ��
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

% �����ǵװ� label_and_eͼƬ
fig_path = [ trackpath, '\GT\label_and_e\'];
fig_dir = dir([ fig_path, '*.fig' ]);  
frame = numel(fig_dir);

% ����ֱ�Ӽ���õ�label2e�����û������¼���
label2e_path = [ trackpath, '\GT\Label_to_Ellipse.mat'];
if exist(label2e_path,'file')
    load(label2e_path, 'label2e');
else
    generate_label2e(); % ����label2eͼƬ
end

% ����superpixel
load([trackpath, '\Pair\Pre_data_New.mat'],'SuperPixel');

% ����ȫ�ֱ��� GT_move �洢��׼Ǩ�ƴ�
% ����ȫ�ֱ��� GT_delete �洢��ǰ�������ҪҪ��ɾȥ���¼�
global GT_move GT_delete global_x global_y t;

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

%% ���б��

for t = 2:2 %frame-1 % 2015.12.14��fluoѵ����������
    
    % �������2������ֹ�����޸�ʱ����
    GT_move{t} = cell(100,5); % Ԥ��20��Ӧ�ù��õ�
    GT_delete{t} = {};
    
    % ��2������ͼƬ
    for ss = t:t+1
        open_fig_and_color_not_labeled( label2e, SuperPixel, fig_dir, ss );
    end
    % ȡ��figure����Բ�ľ��e1,e2
    h1 = get(1, 'children');
    h2 = get(2, 'children');
    e1 = findobj(h1, 'Type', 'Line');
    e2 = findobj(h2, 'Type', 'Line');

    % ���ð�����Ӧ�͵����Ӧ
    global_x = 1;
    global_y = 1;
    set(1, 'keypressfcn', @press_c_space); % ���ո��������
    set(2, 'keypressfcn', @press_c_space); % ��cȷ�����csp��baspѡ�����
    set(e1, 'ButtonDownFcn', @myclick);
    set(e2, 'ButtonDownFcn', @myclick);

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






















