function [ cost_for_train PRF COUNT ] = CX_Calculate_Loss( dataset, exist_GT, s_frame, e_frame, fij, fit, fid, fiv, fmj, fsj )

% 计算精度
%% 读入数据

[ ~, trackpath ] = getpath( dataset );

load([ trackpath, '\Pair\Pre_data_New.mat'], 'n','conflict_pair_last_xy','conflict_pair_next_xy','conflict_fij');
gt_flow_addr = [ trackpath, '\GT\GT_after_hand_tune\GT_Flow_Variables_New.mat'];

% 载入标准答案GT（必须要在GT存在的情况下，使用bool exist_GT来标识GT是否存在）
if exist_GT
    load( gt_flow_addr );
else
    cost_for_train = [];
    PRF = [];
    COUNT = [];
    return;
end

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
    Pcount(1) = Pcount(1) + sum(sum( fij{t})); % 测试中出现的个数
    
    TP.fit = TP.fit + sum(sum( Fit{t}.*fit{t} ));
    Tcount(2) = Tcount(2) + sum(sum( Fit{t}));
    Pcount(2) = Pcount(2) + sum(sum( fit{t}));
    
    TP.fid = TP.fid + sum(sum( Fid{t}.*fid{t} ));
    Tcount(3) = Tcount(3) + sum(sum( Fid{t}));
    Pcount(3) = Pcount(3) + sum(sum( fid{t}));
    
    TP.fiv = TP.fiv + sum(sum( Fiv{t}.*fiv{t} ));
    Tcount(4) = Tcount(4) + sum(sum( Fiv{t}));
    Pcount(4) = Pcount(4) + sum(sum( fiv{t}));
end
for t = s_frame+1:e_frame
    TP.fmj = TP.fmj + sum(sum( Fmj{t}.*fmj{t} ));
    Tcount(5) = Tcount(5) + sum(sum( Fmj{t}));
    Pcount(5) = Pcount(5) + sum(sum( fmj{t}));
    
    TP.fsj = TP.fsj + sum(sum( Fsj{t}.*fsj{t} ));
    Tcount(6) = Tcount(6) + sum(sum( Fsj{t}));
    Pcount(6) = Pcount(6) + sum(sum( fsj{t}));
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
        if 1-all_fin == 1 % 计算测试中虚景出现的次数
            fd_Pcount = fd_Pcount + 1;
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
    if 1-all_fout == 1 % 计算测试中虚景出现的次数
        fd_Pcount = fd_Pcount + 1;
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
cost_for_train = (event_FN + fd_FN)/ sum(Tcount); % 类似recall的计算（仅用于SSVM训练时的cost）

% ------------------------ 下面用于计算测试集精度 ------------------------- %
Preci.all = (event_TP+fd_TP)/(sum(Pcount)+ fd_Pcount);
Preci.move = TP.fij/Pcount(1); % precision=TP/(TP+FP);
Preci.disappear = TP.fit/Pcount(2);
Preci.divide = TP.fid/Pcount(3);
Preci.split = TP.fiv/Pcount(4);
Preci.merge = TP.fmj/Pcount(5);
Preci.appear = TP.fsj/Pcount(6);
Preci.false_detection = fd_TP/fd_Pcount;

Recall.all = (event_TP+fd_TP)/(sum(Tcount)+ fd_Tcount);
Recall.move = (TP.fij)/Tcount(1); % recall=TP/(TP+FN)
Recall.disappear = (TP.fit)/Tcount(2);
Recall.divide = (TP.fid)/Tcount(3);
Recall.split = (TP.fiv)/Tcount(4);
Recall.merge = (TP.fmj)/Tcount(5);
Recall.appear = (TP.fsj)/Tcount(6);
Recall.false_detection = fd_TP/fd_Tcount;

FM.all = 1/(1/Preci.all + 1/Recall.all)*2;
FM.move = 1/(1/Preci.move + 1/Recall.move)*2;
FM.disappear = 1/(1/Preci.disappear + 1/Recall.disappear)*2;
FM.divide = 1/(1/Preci.divide + 1/Recall.divide)*2;
FM.split = 1/(1/Preci.split + 1/Recall.split)*2;
FM.merge = 1/(1/Preci.merge + 1/Recall.merge)*2;
FM.appear = 1/(1/Preci.appear + 1/Recall.appear)*2;
FM.false_detection = 1/(1/Preci.false_detection + 1/Recall.false_detection)*2;

PRF.Preci = Preci; % 用PRF来放置所有的精度相关量
PRF.Recall = Recall;
PRF.FM = FM;

COUNT.Tcount = Tcount; % 用COUNT来放置所有的计数相关量
COUNT.Pcount = Pcount;
COUNT.fd_Tcount = fd_Tcount;
COUNT.fd_Pcount = fd_Pcount;

% ----------------------------------------------------------------------- %


























