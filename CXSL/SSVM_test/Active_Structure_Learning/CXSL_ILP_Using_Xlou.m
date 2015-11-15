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

if 1
    if 1
        disp('  载入之前 SSVM 训练得到的最佳 w...');
        load([ traintrackpath, '\结构化学习\SSVM_Best_W_New.mat']);
    else % 想要复现 excel 中记载的以前的实验结果，只需要手动填写 w_best 即可
        disp('  载入手动填写的w...');
        % 注意！只能用2帧长度样本的训练结果来测试这个！
        
    end
else  
    disp('  载入单独 SVM 训练出来的各事件w，再进行组合...');
    % 也可以使用单独svm训练出来的各事件w，在进行组合
    load([ traintrackpath, '\结构化学习\initial_w_New.mat']);
    % 注意 w 的顺序不能乱
    w_best = [ wij, wit, wid, wiv, wmj, wsj ]';
end

load([ trackpath, '\Pair\Pre_data_New.mat'], 'n');
frame = numel(n);
Fij = cell(frame-1,1);
Fid = cell(frame-1,1);
Fiv = cell(frame-1,1);
Fit = cell(frame-1,1);
Fsj = cell(frame,1);
Fmj = cell(frame,1);

% 调用cplex进行22求解
for i=1:frame-1
    s_frame = i;
    e_frame = i + 1;
    % 使用xlou的 sl 的方法进行训练
    [ Fij, Fit, Fid, Fiv, Fsj, Fmj ] = ASLearning( w_best, s_frame, e_frame, Fij, Fit, Fid, Fiv, Fsj, Fmj );
end


% 调用函数计算精度（假说精度）
s_frame = 1;
e_frame = frame;
[ ~, PRF, COUNT ] = CX_Calculate_Loss( dataset, exist_GT, s_frame, e_frame, Fij, Fit, Fid, Fiv, Fmj, Fsj );

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


% 保存得到的流量变量结果到 track_data 中，供后续画图用（通常不要保存！）
if 0
    txtpath = [ trackpath, '\测试结果记录\BCFWavg_my\40_2_cons1235.txt'];
    fidin = fopen(txtpath, 'w'); fclose(fidin);
    save(strrep(txtpath,'txt','mat'), 'PRF','COUNT','Fij','Fid','Fiv','Fit','Fsj','Fmj'); % 注意修改mat名称
    
%     save([ trackpath, '\结构化学习\Tracking_Data.mat'], 'Fij','Fid','Fiv','Fit','Fsj','Fmj');
end










