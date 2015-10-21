function K = ssvm_kernel_train(psi, alpha, info, kernel_type, cmd,...
    fij, fit, fid, fiv, fmj, fsj,...
    feature_fij_p, feature_fit_p, feature_fid_p, feature_fiv_p, feature_fmj_p, feature_fsj_p)


% 在ssvm训练时，phi1是常数（向量），phi2是变量
% 需要将phi2拆成feature*y，再使用核函数（注意feature是一个矩阵）

switch info.ev % 根据事件选择feature和流程变量
    case 1
        feature = feature_fij_p(info.s_frame:info.e_frame-1);
        flowvar = fij{info.ind}(info.s_frame:info.e_frame-1);
    case 2
        feature = feature_fit_p(info.s_frame:info.e_frame-1);
        flowvar = fit{info.ind}(info.s_frame:info.e_frame-1);
    case 3
        feature = feature_fid_p(info.s_frame:info.e_frame-1);
        flowvar = fid{info.ind}(info.s_frame:info.e_frame-1);
    case 4
        feature = feature_fiv_p(info.s_frame:info.e_frame-1);
        flowvar = fiv{info.ind}(info.s_frame:info.e_frame-1);
    case 5
        feature = feature_fmj_p(info.s_frame:info.e_frame);
        flowvar = fmj{info.ind}(info.s_frame:info.e_frame);
    case 6
        feature = feature_fsj_p(info.s_frame:info.e_frame);
        flowvar = fsj{info.ind}(info.s_frame:info.e_frame);
end

% dims = numel(feature{2}{1});
n_A = size(psi,2);
K = 0; % 相当于线性时<phi1, phi2>的值

for i=1:numel(feature)
    if i==1 && (info.ev==5 || info.ev==6)
        alphaK_mat = [];
    else
        K_mat = zeros(n_A, numel(feature{i})); % 存放psi与feature的核函数的结果
        % psi的每一列与feature的每一列做核运算
        for jj=1:n_A
            tmp = cellfun(@(x)svm_kernel(psi(:,jj), x, kernel_type, cmd), feature{i},'un',0);
            % K_mat为 n_A*n_feature 矩阵
            K_mat(jj,:) = reshape(cell2mat(tmp), 1, []); 
        end
        % 存放aplha筛选完K_mat以后的结果 
        alphaK_mat  = reshape(alpha'*K_mat, size(flowvar{i}));
    end
    K = K + sum(sum(alphaK_mat.*flowvar{i}));
end

% % merge和appear事件第一帧不会发生，要去掉
% if info.ev==5 || info.ev==6
%     K = K - sum(sum(K_cell{1}.*flowvar{1}));
% end % 其他事件最后一轮不会发生（因此没有分配变量空间）此处无需额外去除
    


    
































