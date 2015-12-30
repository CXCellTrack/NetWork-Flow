% ======================================================================= %
%
% 这个是测试 active structured learning 中方法精度的脚本
% 按照2帧长度样本的训练结果w来22测试，得到分配方案F，共需测试frame-1次
% 2015.10.6
% ======================================================================= %

clear;close all

% 指定是否存在GT，如果存在，则计算精度等指标，否则不需要算
exist_GT = 1;
dataset = 'competition';
[ ~, trackpath ] = getpath( dataset );
[ ~, traintrackpath ] = getpath( 'training' );

%% 先求出22之间的约束
global Ellipse n conflict_fij conflict_pair_next_xy conflict_pair_last_xy
load([ trackpath, '\Pair\Pre_data_New.mat']);
frame = numel(n);

fij = cell(frame-1,1);
fid = cell(frame-1,1);
fiv = cell(frame-1,1);
fit = cell(frame-1,1);
fsj = cell(frame,1);
fmj = cell(frame,1);
F = {};

% 预计算约束条件
for t=1:frame-1
    s_frame = t;
    e_frame = t + 1;
    % 使用xlou的 sl 的方法进行训练
    [ fij{t} fit{t} fid{t} fiv{t} fmj{t} fsj{t} F{t} ] = ASLearning( 'ASL', s_frame, e_frame );
end

%% 带入w计算目标函数，并求解（这样分开做易于更换w）
Fij = cell(frame-1,1);
Fid = cell(frame-1,1);
Fiv = cell(frame-1,1);
Fit = cell(frame-1,1);
Fsj = cell(frame,1);
Fmj = cell(frame,1);

% 想要复现 excel 中记载的以前的实验结果，只需要手动填写 w_best 即可
disp('  载入手动填写的w...');
% 注意！只能用2帧长度样本的训练结果来测试这个！  
thisfile = 'FW\ASL\loss_initwp_line_b.mat';
load([ traintrackpath, '\训练结果记录\', thisfile ], 'w_best','Wavg');

if 0
load([ traintrackpath, '\结构化学习\initial_w_New.mat']);
% 注意 w 的顺序不能乱
w_best = [ wij,bij, wit,bit, wid,bid, wiv,biv, wmj,bmj, wsj,bsj ]'; % 原版b
% ------------- 椭圆假说 ------------- %
% 第一个数据集最初w
% w_best = [-24.05027785	0.965959752	-0.209700235	0.023655542	-0.901444678	0.915527485	-0.723055368	0.78127324	-22.34659216	-3.45491283	-1.682414322	-5.355960441	-2.391659001	2.862181421	-7.382944338	8.382838223	1.94377663	-0.451290137	-1.07738777	-4.844423375	-1.122913059	-0.801496889	3.907101647	-11.61160994	3.710115534	0.998335816	4.252699702	0.790594494	1.207125853	3.799458373	1.390618031	5.18991389	1.129864864	0.673380786	-2.076937813	-1.97433464	-1.980221778	-0.051210814	0.597328997	-3.897482158]';
% 第二个数据集最初w
% w_best = [-18.21315239	2.048551055	0.090611096	-0.07830978	-0.681768441	0.091705287	-0.284558766	0.113666465	-16.77170209	-2.820207584	-1.606735489	-2.170929556	-2.000511632	2.42450433	-4.406444861	10.34098417	1.758814312	-0.819906672	-2.159095585	-4.969572233	1.29607646	0.202318113	3.379177651	-14.25716614	7.109097539	3.366559674	2.10659084	-2.499899814	-3.9849466	2.501397849	1.351917771	4.699879458	-0.525793863	-0.261628668	-4.945586679	-1.846739083	-4.998199696	-0.003648138	1.737582541	-8.35761154]';
% 第三个数据集最初w
w_best = [-17.69956928	1.61055833	-0.001787913	0.199592206	-0.287005286	-1.058810184	-0.194715268	0.170821175	-15.86943474	-2.459629861	-1.102503407	-1.233255838	-1.209682446	-2.849135318	-5.932152347	3.635597276	3.112043952	0.905772748	-0.112320949	-0.658275898	4.373948249	0.504445114	-0.211398811	-3.571706443	0	0	0	0	0	0	0	0	0	0	-2.873279295	-0.074622812	-1.438747225	-0.381464188	-2.838608426	-5.896985464]';
% ------------------------------------ %
end

options = sdpsettings('verbose',0,'solver','gurobi');
for t=1:frame-1
    s_frame = t;
    e_frame = t + 1;
    object_function = CXSL_Calculate_Obj_New( dataset, w_best, s_frame, e_frame, fij{t}, fit{t}, fid{t}, fiv{t}, fmj{t}, fsj{t} );
    % 求解
    disp(['  求解 ',num2str(s_frame), '—',num2str(e_frame), ' 帧的目标函数...']);
    sol = solvesdp( F{t}, -object_function, options );
    if sol.problem == 0
        for i = s_frame:e_frame-1
            Fij{i} = round(value(fij{t}{i})) ;
            Fid{i} = round(value(fid{t}{i})) ;
            Fiv{i} = round(value(fiv{t}{i})) ;
            Fit{i} = round(value(fit{t}{i})) ;
        end
        for i = s_frame+1:e_frame
            Fsj{i} = round(value(fsj{t}{i})) ;
            Fmj{i} = round(value(fmj{t}{i})) ;
        end

        COST = value(object_function);
        fprintf('\tcost:\t%.4f\n\n', COST);
        % ------------------------------------------------------ %
    else
        sol.info
        yalmiperror(sol.problem)
    end
end

%% 调用函数计算精度（假说精度）
s_frame = 1;
e_frame = frame;
addfd = 0; % 不在总精度中加入虚景
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
        matpath = [trackpath, '\测试结果记录\local.mat'];
    end
    save(matpath, 'PRF','COUNT','Fij','Fit','Fid','Fiv','Fmj','Fsj'); % 注意修改mat名称
    file = fopen(strrep(matpath,'mat','txt'), 'w'); fclose(file);
end










