function [ fij fit fid fiv fmj fsj ] = CXSL_Assign_FlowVar( dataset, s_frame, e_frame )

% ======================================================================= %
% 这个函数用于分配各事件流程变量，并计算损失函数 sum_cost
% 被 CXSL_ILP_Using_Best_W 直接调用所调用

% 在给定 w 的情况下，生成变量矩阵和目标函数
% 利用公式（2)：L(x,z,w) = w'*fai(x,z)
% 定义 w' = [ wij, wit, wid, wiv, wmj, wsj ]  
% fai(x,z)也按相同的顺序排列
% 最后输出各个事件的变量矩阵和目标函数 object_f，送到外部进行ILP求解
%
% s_frame 为开始帧  e_frame 为结束帧

% ======================================================================= %
% tic;
[ ~, trackpath ] = getpath( dataset );
load([ trackpath, '\Pair\Pre_data_New.mat'], 'n');

%% 2.构建变量矩阵
fij = cell(e_frame-1,1);
fid = cell(e_frame-1,1);
fit = cell(e_frame-1,1);
fiv = cell(e_frame-1,1);
fsj = cell(e_frame,1);
fmj = cell(e_frame,1);


for t = s_frame:e_frame-1
    %  t中第m个前景中的第width个细胞   
    fij{t} = binvar(n(t), 4, 'full'); %%fij变量矩阵
    fit{t} = binvar(n(t), 1, 'full');   %%消失
    fid{t} = binvar(n(t), 6, 'full');   %%母细胞
    fiv{t} = binvar(n(t), 6, 'full');   %%分裂
end

for t = s_frame+1:e_frame
    fmj{t} = binvar(n(t), 6,  'full');   %%融合 ##注意这个矩阵代表（t * t-1），其他都是（t * t+1）
    fsj{t} = binvar(n(t), 1, 'full');   %%出现 
end

% 目标函数 object_f
% w = ones(42,1);
% ============================================
% 新改动：使用这个函数只计算 fai_x_z 和 sum_cost，组合的过程放到主函数里
% object_f = w'* fai_x_z + sum_cost;
% toc;

end

function sum_f = cell_dot_mutil( feature, z )
%
% f 为特征向量cell	n(t)*n(t+1)
% g 为流程变量矩阵
% 输出 sum_f 为该cell内所有特征向量之和
% ===================================

% 1. for 循环式计算 4.75秒
% sum_f = zeros( size(feature{1,1}) );
% [h w] = size(feature);
% for i=1:h
%     for j=1:w
%         sum_f = sum_f + feature{i, j}* z(i,j);
%     end
% end

% 2. arrayfun运算 3.10秒
sum_f = 0;
ss = numel(feature);
% feature1 = cellfun(@(x)x', feature, 'un',0);
feature2 = reshape(feature, ss, 1);
z1 = reshape(z, ss, 1);
f_z = arrayfun(@(x)feature2{x}*z1(x), 1:ss, 'un',0);
for x=1:ss
    sum_f = sum_f + f_z{x};
end

end











