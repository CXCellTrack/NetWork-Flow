function K = ssvm_kernel_train(psi, alpha, info, kernel_type, cmd,...
    fij, fit, fid, fiv, fmj, fsj,...
    feature_fij_p, feature_fit_p, feature_fid_p, feature_fiv_p, feature_fmj_p, feature_fsj_p)


% ��ssvmѵ��ʱ��phi1�ǳ�������������phi2�Ǳ���
% ��Ҫ��phi2���feature*y����ʹ�ú˺�����ע��feature��һ������

switch info.ev % �����¼�ѡ��feature�����̱���
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
K = 0; % �൱������ʱ<phi1, phi2>��ֵ

for i=1:numel(feature)
    if i==1 && (info.ev==5 || info.ev==6)
        alphaK_mat = [];
    else
        K_mat = zeros(n_A, numel(feature{i})); % ���psi��feature�ĺ˺����Ľ��
        % psi��ÿһ����feature��ÿһ����������
        for jj=1:n_A
            tmp = cellfun(@(x)svm_kernel(psi(:,jj), x, kernel_type, cmd), feature{i},'un',0);
            % K_matΪ n_A*n_feature ����
            K_mat(jj,:) = reshape(cell2mat(tmp), 1, []); 
        end
        % ���aplhaɸѡ��K_mat�Ժ�Ľ�� 
        alphaK_mat  = reshape(alpha'*K_mat, size(flowvar{i}));
    end
    K = K + sum(sum(alphaK_mat.*flowvar{i}));
end

% % merge��appear�¼���һ֡���ᷢ����Ҫȥ��
% if info.ev==5 || info.ev==6
%     K = K - sum(sum(K_cell{1}.*flowvar{1}));
% end % �����¼����һ�ֲ��ᷢ�������û�з�������ռ䣩�˴��������ȥ��
    


    
































