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
%     [ Fij, Fit, Fid, Fiv, Fsj, Fmj ] = ASLearning( w_best, s_frame, e_frame, Fij, Fit, Fid, Fiv, Fsj, Fmj );
    [ fij{t} fit{t} fid{t} fiv{t} fmj{t} fsj{t} F{t} ] = ASLearning( s_frame, e_frame );
end

%% 带入w计算目标函数，并求解（这样分开做易于更换w）
Fij = cell(frame-1,1);
Fid = cell(frame-1,1);
Fiv = cell(frame-1,1);
Fit = cell(frame-1,1);
Fsj = cell(frame,1);
Fmj = cell(frame,1);

if 0
    disp('  载入之前 SSVM 训练得到的最佳 w...');
    load([ traintrackpath, '\结构化学习\SSVM_Best_W_New.mat']);
else % 想要复现 excel 中记载的以前的实验结果，只需要手动填写 w_best 即可
    disp('  载入手动填写的w...');
    % 注意！只能用2帧长度样本的训练结果来测试这个！  
    thisfile = 'BCFWavg_paper\40_2\loss_40_2_cons35_cost1_initwp_line.mat';
    load([ traintrackpath, '\训练结果记录\', thisfile ], 'w_best','Wavg');
%     w_best = Wavg{1792};
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










