function KKK = ssvm_kernel_train_paper(n_SV, alpha_i, info, kernel_type, cmd,...
                    fij, fit, fid, fiv, fmj, fsj, s_frame, e_frame,...
                    feature_fij, feature_fit, feature_fid, feature_fiv, feature_fmj, feature_fsj)    

load([ trackpath, '\GT\GT_after_hand_tune\GT_Flow_Variables_New.mat']); % 载入标准答案
% 在ssvm训练时，phi1是常数（向量），phi2是变量
% 需要将phi2拆成feature*y，再使用核函数（注意feature是一个矩阵）
y_star_feature = cell(info.N,1); % y*对应的特征
y_star = cell(info.N,1); % y*

switch info.ev % 根据事件选择feature和流程变量
    case 1
        for ind=1:info.N
            y_star_feature{ind} = feature_fij(s_frame(ind):e_frame(ind)-1);
            y_star{ind} = Fij(s_frame(ind):e_frame(ind)-1);
        end
        sample_feature = feature_fij(info.s_frame:info.e_frame-1);
        flowvar = fij{info.ind}(info.s_frame:info.e_frame-1);
    case 2
        for ind=1:info.N
            y_star_feature{ind} = feature_fit(s_frame(ind):e_frame(ind)-1);
            y_star{ind} = Fit(s_frame(ind):e_frame(ind)-1);
        end
        sample_feature = feature_fit(info.s_frame:info.e_frame-1);
        flowvar = fit{info.ind}(info.s_frame:info.e_frame-1);
    case 3
        for ind=1:info.N
            y_star_feature{ind} = feature_fid(s_frame(ind):e_frame(ind)-1);
            y_star{ind} = Fid(s_frame(ind):e_frame(ind)-1);
        end
        sample_feature = feature_fid(info.s_frame:info.e_frame-1);
        flowvar = fid{info.ind}(info.s_frame:info.e_frame-1);
    case 4
        for ind=1:info.N
            y_star_feature{ind} = feature_fiv(s_frame(ind):e_frame(ind)-1);
            y_star{ind} = Fiv(s_frame(ind):e_frame(ind)-1);
        end
        sample_feature = feature_fiv(info.s_frame:info.e_frame-1);
        flowvar = fiv{info.ind}(info.s_frame:info.e_frame-1);
    case 5
        for ind=1:info.N
            y_star_feature{ind} = feature_fmj(s_frame(ind)+1:e_frame(ind));
            y_star{ind} = Fmj(s_frame(ind)+1:e_frame(ind));
        end
        sample_feature = feature_fmj(info.s_frame+1:info.e_frame);
        flowvar = fmj{info.ind}(info.s_frame+1:info.e_frame);
    case 6
        for ind=1:info.N
            y_star_feature{ind} = feature_fsj(s_frame(ind)+1:e_frame(ind));
            y_star{ind} = Fsj(s_frame(ind)+1:e_frame(ind));
        end
        sample_feature = feature_fsj(info.s_frame+1:info.e_frame);
        flowvar = fsj{info.ind}(info.s_frame+1:info.e_frame);
end



KKK = 0;
tic;
for ind=1:info.N
    disp(['  计算样本',num2str(ind), '...'])
    % fe1为y*特征，fe2为样本特征
    fe1 = y_star_feature{ind};
    y1 = y_star{ind};
    fe2 = sample_feature;
    y2 = flowvar;
    
    frame_result = binvar(numel(fe1),numel(fe2));
    for ii=1:numel(fe1)
        for jj=1:numel(fe2) % ii,jj表示是前者第i帧和后者第j帧
            
            K_mat = zeros(numel(fe1{ii}),numel(fe2{jj}));
            for kk=1:numel(fe1{ii})
                for mm=1:numel(fe2{jj})
                    K_mat(kk,mm) = svm_kernel(fe1{ii}{kk}, fe2{jj}{mm}, 'linear', ''); %kernel_type, cmd
                end
            end
            % 将y1和y2拉直
            tmpy1 = reshape(y1{ii}, 1, []);
            tmpy2 = reshape(y2{jj}, [], 1);
            % 这个相当于 ΣiΣj( yi*yj*k(fei,fej) )
            frame_result(ii,jj) = tmpy1*K_mat*tmpy2;
        end
    end
    
    % 这个和还要*该样本中支持向量的个数
    KKK = KKK + n_SV(ind, info.ev)*sum(frame_result(:));

end
toc
% 这里面包含了5重循环：样本、前者帧、后者帧、前者向量、后者向量，计算非常耗时！
% 求<f*, f>无需*alpha，因为alpha和一定为1













% dims = numel(feature{2}{1});
n_A = size(psi,2);
K = 0; % 相当于线性时<phi1, phi2>的值

for i=1:numel(sample_feature)
    if i==1 && (info.ev==5 || info.ev==6)
        alphaK_mat = [];
    else
        K_mat = zeros(n_A, numel(sample_feature{i})); % 存放psi与feature的核函数的结果
        % psi的每一列与feature的每一列做核运算
        for jj=1:n_A
            tmp = cellfun(@(x)svm_kernel(psi(:,jj), x, kernel_type, cmd), sample_feature{i},'un',0);
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
    


    
































