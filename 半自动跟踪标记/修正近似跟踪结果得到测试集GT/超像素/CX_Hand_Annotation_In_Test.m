% ================================================================== %
%
% CX 2015.12.13
% 这个脚本用于从原始的超像素假说中标记出正确答案
%
% 标记依据

% 主要步骤包括：1 无中生有类型的，需要删去错误轨迹
%              2 可以修改的，则修改为新的正确轨迹
% GUI使用方法：
%
%   1 迁移：第一帧中点击几个bsp组成csp，第二帧中点击几个bsp组成csp，按空格键确定
%   2 消失：第一帧中点击几个bsp组成csp，按空格键确定
%   3 出现：第二帧中点击几个bsp组成csp，按空格键确定
%   4 分裂：第一帧中点击几个bsp组成csp母细胞，第二帧选完第一个子细胞需要按c确定，再选第二个，最后按空格键确定
%   5 合并：第一帧选完第一个源细胞需要按c确定，再选第二个，第二帧中点击几个bsp组成csp大细胞，最后按空格键确定
%
%    Attention：如果出现其他形式的操作，将会报错！
% ================================================================== %
clear;close all;
if 0
    dataset = 'competition';
else
    dataset = 'training';
end
[ ~, trackpath ] = getpath( dataset );

fig_path = [ trackpath, '\???\'];
if strcmp(dataset, 'training')
    fig_path = [ trackpath, '\GT\label_and_e\'];
end
screen_size = get(0,'ScreenSize');
fig_dir = dir([ fig_path, '*.fig' ]);  
frame = numel(fig_dir);


% 定义全局变量 GT_move 存储标准迁移答案
% 定义全局变量 GT_delete 存储当前标记中需要要被删去的事件
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


for t = 1:5 %frame-1 %

    % 先清空这2个，防止后续修改时出错
    GT_move{t} = cell(100,4); % 预留20行应该够用的
    GT_delete{t} = {};
    
    % 打开2副相邻图片
    for ss = t:t+1
        fig_name = [ fig_path, fig_dir(ss).name ];
        openfig(fig_name, 'new', 'visible');
        set(gcf, 'Position', screen_size);
    end
    % 取得figure内椭圆的句柄e1,e2
    h1 = get(1, 'children');
    h2 = get(2, 'children');
    e1 = findobj(h1, 'Type', 'Line');
    e2 = findobj(h2, 'Type', 'Line');

    % 设置按键响应和点击响应
    global_x = 1;
    global_y = 1;
    set(1, 'keypressfcn', @press_c_space); % 按空格键负责换行
    set(2, 'keypressfcn', @press_c_space); % 按c确认组合csp的basp选择完毕
    set(e1, 'ButtonDownFcn', @myclick_1);
    set(e2, 'ButtonDownFcn', @myclick_2);

    disp(['  正在标记',num2str(t),'—',num2str(t+1),'帧']);
    
    % 通过input操作来给标记留出时间
    keyin = input('  在控制台上按enter键进入下两帧\n', 's');
    if strcmp(keyin, '')
        close(1);
        close(2);
    end
    
    if 1
        disp('  保存结果到Hand_GT_New.mat 中');
        GT_delete_s = GT_delete;
        GT_move_s = GT_move;
        save([ trackpath, '\GT\Hand_GT_New.mat'], 'GT_move_s','GT_delete_s');
    end

end






















