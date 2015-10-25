function K = ssvm_kernel_train_my(psi, alpha, info, kernel_type_ev, cmd_ev,...
    fij, fit, fid, fiv, fmj, fsj)

% 载入特征
[ ~, trackpath ] = getpath( 'training' );
load([ trackpath, '\结构化学习\Feature_New.mat']); 
% 在ssvm训练时，phi1是常数（向量），phi2是变量
% 需要将phi2拆成feature*y，再使用核函数（注意feature是一个矩阵）

switch info.ev % 根据事件选择feature和流程变量
    case 1
        feature = feature_fij(info.s_frame:info.e_frame-1); % 训练特征
        flowvar = fij{info.ind}(info.s_frame:info.e_frame-1); % 训练变量
    case 2
        feature = feature_fit(info.s_frame:info.e_frame-1);
        flowvar = fit{info.ind}(info.s_frame:info.e_frame-1);
    case 3
        feature = feature_fid(info.s_frame:info.e_frame-1);
        flowvar = fid{info.ind}(info.s_frame:info.e_frame-1);
    case 4
        feature = feature_fiv(info.s_frame:info.e_frame-1);
        flowvar = fiv{info.ind}(info.s_frame:info.e_frame-1);
    case 5
        feature = feature_fmj(info.s_frame+1:info.e_frame);
        flowvar = fmj{info.ind}(info.s_frame+1:info.e_frame);
    case 6
        feature = feature_fsj(info.s_frame+1:info.e_frame);
        flowvar = fsj{info.ind}(info.s_frame+1:info.e_frame);
end

% 支持向量个数n_A
n_A = size(psi,2);
K = 0; % 相当于线性时<phi1, phi2>的值

for tt=1:numel(feature) % 按帧循环
    
    K_mat = zeros(n_A, numel(feature{tt})); % 存放psi与feature的核函数的结果
    % psi的每一列与feature的每一列做核运算
    for jj=1:n_A
        % 此处使用for循环和使用cellfun的时间是一样的
        tmp = cellfun(@(x)svm_kernel(psi(:,jj), x, kernel_type_ev, cmd_ev), feature{tt},'un',0);
        % K_mat为 n_A*n_feature 矩阵
        K_mat(jj,:) = reshape(cell2mat(tmp), 1, []); 
    end
    % 存放aplha筛选完K_mat以后的结果 
    alphaK_mat  = reshape(alpha'*K_mat, size(flowvar{tt}));
    K = K + sum(sum(alphaK_mat.*flowvar{tt}));
end

    


    
































