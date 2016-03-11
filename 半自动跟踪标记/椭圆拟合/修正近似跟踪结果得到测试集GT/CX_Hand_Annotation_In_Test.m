% ================================================================== %
%
% CX 2015.7.21
% 这个脚本用于将一个较为接近gt的trackdata figure手动修正为ground truth figure
% 较为接近gt的data来自于训练集得到的w的在测试集上的结果
%
% 标记依据

% 主要步骤包括：1 无中生有类型的，需要删去错误轨迹
%              2 可以修改的，则修改为新的正确轨迹
% GUI使用方法：
%
%   1 选中椭圆，按d键，删除此椭圆的来路和去路（只能在前一帧中操作），按空格键确定
%   2 选中第一帧椭圆，再选中第二帧椭圆，添加一条迁移轨迹，按空格键确定
%   3 选中第一帧椭圆，按空格键确定，添加消失事件
%   4 直接选中第二帧椭圆，按空格键确定，添加出现事件
%   5 选中第一帧的1个，再选中第二帧的2个，按空格键确定，添加divide/split事件
%   6 选中第一帧的2个，再选中第二帧的1个，按空格键确定，添加merge事件
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

fig_addr = [ trackpath, '\Pair\可视化跟踪标记\'];
if strcmp(dataset, 'training')
    fig_addr = [ trackpath, '\GT\'];
end
screen_size = get(0,'ScreenSize');
colored_fig_dir = dir([ fig_addr, '*.fig' ]);  
frame = numel(colored_fig_dir);


% 定义全局变量 GT_move 存储标准迁移答案
% 定义全局变量 GT_delete 存储当前标记中需要要被删去的事件
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


for t = 1:frame-1 % 2015.10.4将hela1训练集的前80帧标记完毕

    % 先清空这2个，防止后续修改时出错
    GT_move{t} = zeros(50,4); % 预留20行应该够用的
    GT_delete{t} = [];
    
    % 打开2副相邻图片
    for ss = t:t+1
        fig_name = [ fig_addr, colored_fig_dir(ss).name ];
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
    set(1, 'keypressfcn', @pressBlank); % 按空格键负责换行
    set(2, 'keypressfcn', @pressBlank);
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






















