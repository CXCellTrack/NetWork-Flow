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
% 指定是否存在GT，如果存在，则计算精度等指标，否则不需要算
exist_GT = 1;
disp(['  计算 ',num2str(s_frame), '—',num2str(e_frame), ' 帧的目标函数和约束条件...']);

%% 计算约束条件 注意：不包含损失函数
% ----------------------------------------------------------------------- %
% B方法：先计算 prob = <w,feature>（41s）在计算 obj = prob.*z（1s）
%        如果放在循环中，每次花一秒计算obj略长，但只计算一次的话速度很快
% 因此此处使用B方法速度较快，经验证 count_F_false 的计算没有问题
disp('分配流程变量...');tic
[ fij fit fid fiv fmj fsj ] = CXSL_Assign_FlowVar( dataset, s_frame, e_frame );
toc;disp('计算约束条件...');
% 此处的true/false决定是否加入可选约束（要与训练时的选择一致）
use_op_cons = [3 5];
[ F ] = CXSL_Calculate_Constraint_New_Conflict( dataset, use_op_cons, s_frame, e_frame, fij, fit, fid, fiv, fmj, fsj);
% 计算目标函数（需要载入之前计算好的特征）

%% 组建目标函数
[ ~, traintrackpath ] = getpath( 'training' );

disp('  载入之前 SSVM 训练保存数据...'); 
thisfile = '\核记录\BCFWavg\loss_5_13_init0p_noline-linear.mat';
load([ traintrackpath, '\训练结果记录\', thisfile]);

disp('组建目标函数...');tic
object_function = CXSL_Calculate_Obj_hunhe_paper( dataset, s_frame, e_frame,...
    w_best, delta_y_best, feature_best, alpha_best, kernel_type, cmd, lambda, islinear, ...
    fij, fit, fid, fiv, fmj, fsj );
toc

%% 最终求解
disp('  开始求解ILP...');
clearvars -except F object_function s_frame e_frame  fij fid fiv fit fsj fmj loss dataset count count_F_false exist_GT trackpath thisfile;
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
fprintf('\tcost:\t%.4f\n\n', COST);

%% 调用函数计算精度（假说精度）
addfd = 0; % 选择是否将虚景计算在总精度中
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
if 1
    if exist('thisfile', 'var')
        matpath = [trackpath, '\测试结果记录\',strrep(thisfile, 'loss_', '')];
    else
        matpath = [trackpath, '\测试结果记录\???.mat'];
    end
    save(matpath, 'PRF','COUNT','Fij','Fit','Fid','Fiv','Fmj','Fsj'); % 注意修改mat名称
    file = fopen(strrep(matpath,'mat','txt'), 'w'); fclose(file);
end






