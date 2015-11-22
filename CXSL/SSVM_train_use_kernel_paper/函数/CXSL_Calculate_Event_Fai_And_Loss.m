function [ fij fit fid fiv fmj fsj phi_x_z cost_for_train cost_for_train_all ] = CXSL_Calculate_Event_Fai_And_Loss( s_frame, e_frame )

%
% 在给定 w 的情况下，生成变量矩阵和目标函数
% 利用公式（2)：L(x,z,w) = w'*fai(x,z)
% 定义 w' = [ wij, wit, wid, wiv, wmj, wsj ]  
% fai(x,z)也按相同的顺序排列
% 最后输出各个事件的变量矩阵和目标函数 object_f，送到外部进行ILP求解
%
% s_frame 为开始帧  e_frame 为结束帧
% tic;
dataset = 'training'; % 这个函数也只在训练中使用
[ ~, trackpath ] = getpath( dataset );

load([ trackpath, '\Pair\Pre_data_New.mat'], 'n','conflict_pair_last_xy','conflict_pair_next_xy','conflict_fij');
% 载入特征
load([ trackpath, '\结构化学习\Feature_Plus_New.mat']);
% 载入标准答案GT
load([ trackpath, '\GT\GT_after_hand_tune\GT_Flow_Variables_New.mat']);

% 步骤1已被注释
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

%% 3.生成原始目标函数
% 利用公式（2)：L(x,z,w) = w'*fai(x,z)
% 利用公式（8): delta(z,z*) = 1/|z*| * sum(z* * (1-z))

fai_fij = 0;
fai_fit = 0;
fai_fid = 0;
fai_fiv = 0;
fai_fmj = 0;
fai_fsj = 0;

% ======================================================================= %
% 在特征不全的情况下要使用步骤1中得到映射表计算
% 在特征齐全的情况下直接矩阵点乘更快（6.5秒 VS 31.4秒）
% ======================================================================= %
% A.使用矩阵乘法计算：
% tic;
for t = s_frame:e_frame-1
    % 去掉sum后之间相加速度更快 由 2.9168s 到 2.82s
    fai_fij = fai_fij + cell_dot_mutil( feature_fij_p{t}, fij{t} );

    fai_fit = fai_fit + cell_dot_mutil( feature_fit_p{t}, fit{t} );

    fai_fid = fai_fid + cell_dot_mutil( feature_fid_p{t}, fid{t} );

    fai_fiv = fai_fiv + cell_dot_mutil( feature_fiv_p{t}, fiv{t} );
end
% toc;
for t = s_frame+1:e_frame
    
    fai_fmj = fai_fmj + cell_dot_mutil( feature_fmj_p{t}, fmj{t} );

    fai_fsj = fai_fsj + cell_dot_mutil( feature_fsj_p{t}, fsj{t} );
end

% 将所有的 fai_f 组合成列向量（增广向量）
% fai_x_z = [ fai_fij; fai_fit; fai_fid; fai_fiv; fai_fmj; fai_fsj; ];
% 在使用核时还是用分开事件的 phi 比较灵活
phi_x_z = cell(6,1);
phi_x_z{1} = fai_fij;
phi_x_z{2} = fai_fit;
phi_x_z{3} = fai_fid;
phi_x_z{4} = fai_fiv;
phi_x_z{5} = fai_fmj;
phi_x_z{6} = fai_fsj;

%% 4.计算损失函数
% ======================================================================= %
% 求解各项事件的损失函数，根据公式（10），这部分要放在目标函数里
% ======================================================================= %
% ============== delta(f,f*) ============== %
TP.fij = 0;
TP.fit = 0;
TP.fid = 0;
TP.fiv = 0;
TP.fmj = 0;
TP.fsj = 0;

% PN代表测试结果，TF代表GT结果
%    T  F
% P  TP FP
% N  FN TN
% 
% precision = tp/(tp+fp)
% recall = tp/(tp+fn);

Tcount = zeros(6,1); % tp+tn
Pcount = zeros(6,1); % tp+fn
% 损失函数都按照 false negative 计算，即 f*=1 && f=0 时才有损失
for t = s_frame:e_frame-1
    TP.fij = TP.fij + sum(sum( Fij{t}.*fij{t} ));
    Tcount(1) = Tcount(1) + sum(sum( Fij{t})); % 统计出 Fij{t}中1的个数
    
    TP.fit = TP.fit + sum(sum( Fit{t}.*fit{t} ));
    Tcount(2) = Tcount(2) + sum(sum( Fit{t}));
    
    TP.fid = TP.fid + sum(sum( Fid{t}.*fid{t} ));
    Tcount(3) = Tcount(3) + sum(sum( Fid{t}));
    
    TP.fiv = TP.fiv + sum(sum( Fiv{t}.*fiv{t} ));
    Tcount(4) = Tcount(4) + sum(sum( Fiv{t}));
end
for t = s_frame+1:e_frame
    TP.fmj = TP.fmj + sum(sum( Fmj{t}.*fmj{t} ));
    Tcount(5) = Tcount(5) + sum(sum( Fmj{t}));
    
    TP.fsj = TP.fsj + sum(sum( Fsj{t}.*fsj{t} ));
    Tcount(6) = Tcount(6) + sum(sum( Fsj{t}));
end

% ============= 注意：论文中的方法是吧所有的count加起来在计算 ============== %
event_TP = TP.fij + TP.fit + TP.fid + TP.fiv + TP.fmj + TP.fsj;
 
%% 计算（GT中不采用的假说 测试却被采用）造成的误差
% ======================================================================= %
% 对于虚景或者矛盾假说的情况，即所有 f* 全为0，但f中出现了1，需要额外计算
% ======================================================================= %

% 考虑2到最后帧的入口，若椭圆j真实入口为0，且分配入口为1，则记作一次损失
fd_TP = 0;
fd_Pcount = 0; % 统计测试中“虚景“出现次数
fd_Tcount = 0; % 统计GT中“虚景”（即入口出口为0的椭圆）出现的次数

for t = s_frame+1:e_frame
    for j=1:n(t)  
        % 分配入口变量
        sum_fid = 0;
        sum_fiv = 0;
        % 真实入口变量
        sum_Fid = 0;
        sum_Fiv = 0;
        % sum_fid 为所有分裂到包含 j 的 pair 的 fid 之和
        for ind=1:numel(conflict_pair_next_xy{t}{j})/2
            sum_fid = sum_fid + fid{t-1}( conflict_pair_next_xy{t}{j}(ind,1), conflict_pair_next_xy{t}{j}(ind,2) );
            sum_fiv = sum_fiv + fiv{t-1}( conflict_pair_next_xy{t}{j}(ind,1), conflict_pair_next_xy{t}{j}(ind,2) );
            
            sum_Fid = sum_Fid + Fid{t-1}( conflict_pair_next_xy{t}{j}(ind,1), conflict_pair_next_xy{t}{j}(ind,2) );
            sum_Fiv = sum_Fiv + Fiv{t-1}( conflict_pair_next_xy{t}{j}(ind,1), conflict_pair_next_xy{t}{j}(ind,2) );
        end
        % 该椭圆分配到的入口和 fin　只能为０或１（由约束限制）
        % ------------------------------------- %
        % 用sum_fij代替原来的sum(fij{t-1}(:,j))
        sum_fij = 0;
        % 用sum_Fij代替原来的sum(Fij{t-1}(:,j))
        sum_Fij = 0;
        for ind=1:size(conflict_fij{t-1}{j}, 1)           
            sum_fij = sum_fij + fij{t-1}( conflict_fij{t-1}{j}(ind,1), conflict_fij{t-1}{j}(ind,2) );
            sum_Fij = sum_Fij + Fij{t-1}( conflict_fij{t-1}{j}(ind,1), conflict_fij{t-1}{j}(ind,2) );
        end
        % ------------------------------------- %
        all_fin = sum_fij + fsj{t}(j) + sum(fmj{t}(j,:)) + sum_fid + sum_fiv;        
        % 该椭圆的标准答案 入口和 Fin　（同上）
        all_Fin = sum_Fij + Fsj{t}(j) + sum(Fmj{t}(j,:)) + sum_Fid + sum_Fiv;
        % ------------------------------------- %   
        
        if 1-all_Fin == 1 % 计算GT中虚景出现的次数
            fd_Tcount = fd_Tcount + 1;
        end

        % 定义损失函数，只有当该椭圆真实入口和为0（即虚景），而分配入口和为1时，才算做损失
        fd_TP = fd_TP + (1 - all_Fin)*(1 - all_fin);
        
    end
end

% ===========考虑第一帧的出口，若椭圆j真实出口为0，且分配出口为1，则记作一次损失
t = s_frame;
for j=1:n(t)
    sum_fmj = 0;
    sum_Fmj = 0;
    
    % sum_fmj 为所有包含 j 的融合 pair 的 fmj 之和
    for ind=1:numel(conflict_pair_last_xy{t}{j})/2
        sum_fmj = sum_fmj + fmj{t+1}( conflict_pair_last_xy{t}{j}(ind,1), conflict_pair_last_xy{t}{j}(ind,2) );
        sum_Fmj = sum_Fmj + Fmj{t+1}( conflict_pair_last_xy{t}{j}(ind,1), conflict_pair_last_xy{t}{j}(ind,2) );
    end
    % 该椭圆分配到的出口和 fout　只能为０或１（由约束限制）
    all_fout = sum(fij{t}(j,:)) + fit{t}(j) + sum(fid{t}(j,:)) + sum(fiv{t}(j,:)) + sum_fmj;
    % 该椭圆的标准答案出口和 Fout　（同上）
    all_Fout = sum(Fij{t}(j,:)) + Fit{t}(j) + sum(Fid{t}(j,:)) + sum(Fiv{t}(j,:)) + sum_Fmj;
            
    if 1-all_Fout == 1 % 计算GT中虚景出现的次数
        fd_Tcount = fd_Tcount + 1;
    end

    % 定义损失函数，只有当该椭圆真实出口和为0（即虚景），而分配出口和为1时，才算做损失
    fd_TP = fd_TP + (1 - all_Fout)*(1 - all_fout);
end

% 因为有多假说的存在，“虚景”出现次数必不为0
% 注意，此处所指的“虚景”包括了 真正的虚景 和 被矛盾集排除掉的椭圆，确定规则为：入口/出口是否为0
    
%% 5.得到目标函数
% ======================================================================= %
%
% 总结目标函数的组成：w'* fai(x,z) + delta(z,z*)
% 其中 delta(z,z*) 又由2部分组成：
% (1) F*(1 - f)     文章中的定义    1 0时记作损失 考察每个事件的流程变量f
%     共6个: TN.fij = 0;
%            TN.fit = 0;
%            TN.fid = 0;
%            TN.fiv = 0;
%            TN.fmj = 0;
%            TN.fsj = 0;
%
% (2) (1 - F)* f	我自己的定义   0 1时记作损失 只考察各事件流程变量全0的椭圆
%     共2个: 2到最后帧的入口cost和第一帧的出口cost
%            记作 add_cost
%
% ======================================================================= %
% 全部损失函数cost之和
% 我自己的方法（先注释掉）
% sum_cost = TN.fij + TN.fit + TN.fid + TN.fiv + TN.fmj + TN.fsj + add_cost;
% 论文方法

event_FN = sum(Tcount) - event_TP; % 细胞事件的FN（测试没发生但实际发生了）
fd_FN = fd_Tcount - fd_TP; % 虚景的FN（测试不为虚但实际为虚）
cost_for_train_all = (event_FN + fd_FN)/ sum(Tcount); % 类似recall的计算（仅用于SSVM训练时的cost）
cost_for_train = event_FN/ sum(Tcount);

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
f_z = arrayfun(@(x) feature2{x}*z1(x), 1:ss, 'un',0);
for x=1:ss
    sum_f = sum_f + f_z{x};
end

end











