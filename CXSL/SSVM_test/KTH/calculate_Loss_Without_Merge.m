function [ PRF COUNT ] = calculate_Loss_Without_Merge( dataset, exist_GT, s_frame, e_frame, Pcount, fij, fit, fid, fsj )

% 计算精度
%% 读入数据

[ ~, trackpath ] = getpath( dataset );

load([ trackpath, '\Pair\Pre_data_New.mat'], 'n','conflict_pair_last_xy','conflict_pair_next_xy','conflict_fij');
gt_flow_addr = [ trackpath, '\GT\GT_after_hand_tune\GT_Flow_Variables_New.mat'];

% 载入标准答案GT（必须要在GT存在的情况下，使用bool exist_GT来标识GT是否存在）
if exist_GT
    load( gt_flow_addr );
else
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
TP.fsj = 0;

% PN代表测试结果，TF代表GT结果
%    T  F
% P  TP FP
% N  FN TN
% 
% precision = tp/(tp+fp)
% recall = tp/(tp+fn);

Tcount = zeros(4,1); % tp+tn
% 损失函数都按照 false negative 计算，即 f*=1 && f=0 时才有损失
for t = s_frame:e_frame-1
    TP.fij = TP.fij + sum(sum( Fij{t}.*fij{t} ));
    Tcount(1) = Tcount(1) + sum(sum( Fij{t})); % 统计出 Fij{t}中1的个数
%     Pcount(1) = Pcount(1) + sum(sum( fij{t})); % 测试中出现的个数
    
    TP.fit = TP.fit + sum(sum( Fit{t}.*fit{t} ));
    Tcount(2) = Tcount(2) + sum(sum( Fit{t}));
    
    TP.fid = TP.fid + sum(sum( Fid{t}.*fid{t} ));
    Tcount(3) = Tcount(3) + sum(sum( Fid{t}));
end
for t = s_frame+1:e_frame
    TP.fsj = TP.fsj + sum(sum( Fsj{t}.*fsj{t} ));
    Tcount(4) = Tcount(4) + sum(sum( Fsj{t}));
end

% ============= 注意：论文中的方法是吧所有的count加起来在计算 ============== %
event_TP = TP.fij + TP.fit + TP.fid + TP.fsj;
 
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

% ------------------------ 下面用于计算测试集精度 ------------------------- %
Preci.all = event_TP/(sum(Pcount));
Preci.move = TP.fij/Pcount(1); % precision=TP/(TP+FP);
Preci.disappear = TP.fit/Pcount(2);
Preci.divide = TP.fid/Pcount(3);
Preci.appear = TP.fsj/Pcount(4);

Recall.all = event_TP/(sum(Tcount));
Recall.move = (TP.fij)/Tcount(1); % recall=TP/(TP+FN)
Recall.disappear = (TP.fit)/Tcount(2);
Recall.divide = (TP.fid)/Tcount(3);
Recall.appear = (TP.fsj)/Tcount(4);

FM.all = 1/(1/Preci.all + 1/Recall.all)*2;
FM.move = 1/(1/Preci.move + 1/Recall.move)*2;
FM.disappear = 1/(1/Preci.disappear + 1/Recall.disappear)*2;
FM.divide = 1/(1/Preci.divide + 1/Recall.divide)*2;
FM.appear = 1/(1/Preci.appear + 1/Recall.appear)*2;

PRF.Preci = Preci; % 用PRF来放置所有的精度相关量
PRF.Recall = Recall;
PRF.FM = FM;

COUNT.Tcount = Tcount; % 用COUNT来放置所有的计数相关量
COUNT.Pcount = Pcount;

% ----------------------------------------------------------------------- %


























