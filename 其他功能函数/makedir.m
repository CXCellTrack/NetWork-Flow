% 这个脚本用于制作目录，仅在训练或测试开始前使用一次即可！
clear

%% 训练目录
dataset = 'training';
[ ~, trackpath ] = getpath( dataset );
% iiii = 6;
% trackpath = ['E:\datasets\first_edition\training_datasets\N2DH-SIM\0',num2str(iiii),'_0-00_track'];
path = {};
path{1} = [ trackpath, '\GT'];
path{2} = [ trackpath, '\GT\GT_after_hand_tune'];
path{3} = [ trackpath, '\GT\label_and_e'];
path{4} = [ trackpath, '\Pair'];
path{5} = [ trackpath, '\Pair\可视化跟踪标记'];
path{6} = [ trackpath, '\结构化学习'];
path{7} = [ trackpath, '\新拟合图'];
path{8} = [ trackpath, '\训练结果记录'];

for i=1:numel(path)
    mkdir(path{i});
end

%% 测试目录
dataset = 'competition';
[ ~, trackpath ] = getpath( dataset );
path = {};
path{1} = [ trackpath, '\GT'];
path{2} = [ trackpath, '\GT\GT_after_hand_tune'];
% path{3} = [ trackpath, '\GT\label_and_e'];
path{4} = [ trackpath, '\Pair'];
path{5} = [ trackpath, '\Pair\可视化跟踪标记'];
path{6} = [ trackpath, '\结构化学习'];
path{7} = [ trackpath, '\新拟合图'];
path{8} = [ trackpath, '\测试结果记录'];


for i=1:numel(path)
    mkdir(path{i});
end