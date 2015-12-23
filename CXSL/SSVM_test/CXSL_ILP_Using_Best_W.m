% ================================================================== %
%
% CX 2015.6.23
% 这个脚本用于在给定 w 的情况下寻找测试集合中ILP问题的最优解，即最佳分配方案
% 需要自己指定开始帧 s_frame 和结束帧 e_frame 
%
% ================================================================== %

clear;close all;

%% 指定在哪个数据集上进行计算（training or competition）
if 1
    dataset = 'competition';
    disp('<选定测试集>');
else
    dataset = 'training';
    disp('<选定训练集>');
end
[ segpath, trackpath ] = getpath( dataset );

% 指定测试帧的范围
segdir = dir([ segpath, '\*.tif']);
s_frame = 1;
e_frame = numel(segdir);
disp(['  计算 ',num2str(s_frame), '—',num2str(e_frame), ' 帧的目标函数和约束条件...']);

%% 计算约束条件 注意：不包含损失函数
% ----------------------------------------------------------------------- %
% B方法：先计算 prob = <w,feature>（41s）在计算 obj = prob.*z（1s）
%        如果放在循环中，每次花一秒计算obj略长，但只计算一次的话速度很快
% 因此此处使用B方法速度较快，经验证 count_F_false 的计算没有问题
disp('分配流程变量...');tic
[ fij fit fid fiv fmj fsj ] = CXSL_Assign_FlowVar( dataset, s_frame, e_frame );
toc;disp('计算约束条件...');
% 此处的true/false决定是否加入可选约束（要与训练时的选择一致)
use_op_cons_test = [3 5 4];
[ Ffull, Fbase ] = CXSL_Calculate_Constraint_New_Conflict( dataset, use_op_cons_test, s_frame, e_frame, fij, fit, fid, fiv, fmj, fsj);
% 计算目标函数（需要载入之前计算好的特征）
if 1
    F = Ffull;
else
    F = Fbase; % 不用可选约束
end

%% 组建目标函数
[ ~, traintrackpath ] = getpath( 'training' );
if 1
    if 0
        disp('  载入之前 SSVM 训练得到的最佳 w...');
        load([ traintrackpath, '\结构化学习\SSVM_Best_W_New.mat']);
    else % 想要复现 excel 中记载的以前的实验结果，只需要手动填写 w_best 即可
        disp('  载入手动填写的w...');
        thisfile = 'BCFWavg_paper\withgap\64_2\loss_64_2_cons35_cost1_initwp_line_b_rng.mat';
        load([ traintrackpath, '\训练结果记录\', thisfile ], 'w_best','use_op_cons','Wavg');
%         w_best = Wavg{394};
%         if ~isequal(use_op_cons, use_op_cons_test)
%             error('测试所用的可选约束和训练不一致！');
%         end
    end
else  
    disp('  载入local SVM 训练出来的各事件w...');
    % 也可以使用单独svm训练出来的各事件w，在进行组合
    load([ traintrackpath, '\结构化学习\initial_w_New.mat']);
    % 加偏置
    w_best = [ wij,bij, wit,bit, wid,bid, wiv,biv, wmj,bmj, wsj,bsj ]';
    % 第一个数据集
%     w_best = [-24.05027785	0.965959752	-0.209700235	0.023655542	-0.901444678	0.915527485	-0.723055368	0.78127324	-22.34659216	-3.45491283	-1.682414322	-5.355960441	-2.391659001	2.862181421	-7.382944338	8.382838223	1.94377663	-0.451290137	-1.07738777	-4.844423375	-1.122913059	-0.801496889	3.907101647	-11.61160994	3.710115534	0.998335816	4.252699702	0.790594494	1.207125853	3.799458373	1.390618031	5.18991389	1.129864864	0.673380786	-2.076937813	-1.97433464	-1.980221778	-0.051210814	0.597328997	-3.897482158]';
    % 第二个数据集
%     w_best = [-18.21315239	2.048551055	0.090611096	-0.07830978	-0.681768441	0.091705287	-0.284558766	0.113666465	-16.77170209	-2.820207584	-1.606735489	-2.170929556	-2.000511632	2.42450433	-4.406444861	10.34098417	1.758814312	-0.819906672	-2.159095585	-4.969572233	1.29607646	0.202318113	3.379177651	-14.25716614	7.109097539	3.366559674	2.10659084	-2.499899814	-3.9849466	2.501397849	1.351917771	4.699879458	-0.525793863	-0.261628668	-4.945586679	-1.846739083	-4.998199696	-0.003648138	1.737582541	-8.35761154]';
    % 第三个数据集
%     w_best = [-17.69956928	1.61055833	-0.001787913	0.199592206	-0.287005286	-1.058810184	-0.194715268	0.170821175	-15.86943474	-2.459629861	-1.102503407	-1.233255838	-1.209682446	-2.849135318	-5.932152347	3.635597276	3.112043952	0.905772748	-0.112320949	-0.658275898	4.373948249	0.504445114	-0.211398811	-3.571706443	0	0	0	0	0	0	0	0	0	0	-2.873279295	-0.074622812	-1.438747225	-0.381464188	-2.838608426	-5.896985464]';
end

disp('组建目标函数...');tic
object_function = CXSL_Calculate_Obj_New( dataset, w_best, s_frame, e_frame, fij, fit, fid, fiv, fmj, fsj );
toc

%% 最终求解
disp('  开始求解ILP...');
% clearvars -except F object_function s_frame e_frame  fij fid fiv fit fsj fmj loss dataset count count_F_false exist_GT trackpath thisfile use_op_cons_test;
% 注意，原先采用先算出 fai(x,z) = <feature,z>，在计算 obj = <w,fai(x,z)>;
% 现在采用先计算 <w,feature>，在计算 obj = <w,feature>*z，速度得到了明显提升
% 但这是针对于一次计算而言，如果在循环中每次都要这么计算目标函数，速度还是没有原方法快 2015.6.24

options = sdpsettings('verbose',0,'solver','gurobi');
sol = solvesdp( F, -object_function, options )

Fij = cell(e_frame-1,1);
Fid = cell(e_frame-1,1);
Fiv = cell(e_frame-1,1);
Fit = cell(e_frame-1,1);
Fsj = cell(e_frame,1);
Fmj = cell(e_frame,1);
if sol.problem ~= 0
    sol.info
    yalmiperror(sol.problem);
end

for t = s_frame:e_frame-1
    Fij{t} = round(value(fij{t})) ;
    Fid{t} = round(value(fid{t})) ;
    Fiv{t} = round(value(fiv{t})) ;
    Fit{t} = round(value(fit{t})) ;
end
for t = s_frame+1:e_frame
    Fsj{t} = round(value(fsj{t})) ;
    Fmj{t} = round(value(fmj{t})) ;
end

COST = value(object_function);
fprintf('\tcost:\t%.4f\n\n', COST); % save([trackpath, '\结构化学习\Tracking_Data.mat'], 'Fij','Fid','Fiv','Fit','Fmj','Fsj');

%% 调用函数计算精度（假说精度）
addfd = 0; % 选择是否将虚景计算在总精度中
% 指定是否存在GT，如果存在，则计算精度等指标，否则不需要算
exist_GT = 1;
[ ~, PRF, COUNT ] = CX_Calculate_Loss( dataset, addfd, exist_GT, s_frame, e_frame, Fij, Fit, Fid, Fiv, Fmj, Fsj );

if isa(PRF, 'struct')
    
    disp(PRF.FM) % 只显示F-measure就行了   
    % 打印各事件精度
    Preci = struct2array(PRF.Preci)*100;
    Recall = struct2array(PRF.Recall)*100;
    FMeasure = struct2array(PRF.FM)*100;
    
    PRM_for_excel = [ Preci;Recall;FMeasure ];
    COUNT_for_excel = [ COUNT.Tcount', COUNT.fd_Tcount; COUNT.Pcount', COUNT.fd_Pcount ];
    disp(COUNT_for_excel)
end
% ------------------------------------------------------ %

%% 保存得到的流量变量结果到 track_data 中，供后续画图用（通常不要保存！）
if 0
    if exist('thisfile', 'var')
        matpath = [trackpath, '\测试结果记录\',strrep(thisfile, 'loss_', '')];
    else
        matpath = [trackpath, '\测试结果记录\local_b.mat'];
    end
    save(matpath, 'PRF','COUNT','Fij','Fit','Fid','Fiv','Fmj','Fsj'); % 注意修改mat名称
    file = fopen(strrep(matpath,'mat','txt'), 'w'); fclose(file);
end




