function K = ssvm_kernel_train_my(psi, alpha, info, kernel_type_ev, cmd_ev,...
    fij, fit, fid, fiv, fmj, fsj)

% ��������
[ ~, trackpath ] = getpath( 'training' );
load([ trackpath, '\�ṹ��ѧϰ\Feature_New.mat']); 
% ��ssvmѵ��ʱ��phi1�ǳ�������������phi2�Ǳ���
% ��Ҫ��phi2���feature*y����ʹ�ú˺�����ע��feature��һ������

switch info.ev % �����¼�ѡ��feature�����̱���
    case 1
        feature = feature_fij(info.s_frame:info.e_frame-1); % ѵ������
        flowvar = fij{info.ind}(info.s_frame:info.e_frame-1); % ѵ������
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

% ֧����������n_A
n_A = size(psi,2);
K = 0; % �൱������ʱ<phi1, phi2>��ֵ

for tt=1:numel(feature) % ��֡ѭ��
    
    K_mat = zeros(n_A, numel(feature{tt})); % ���psi��feature�ĺ˺����Ľ��
    % psi��ÿһ����feature��ÿһ����������
    for jj=1:n_A
        % �˴�ʹ��forѭ����ʹ��cellfun��ʱ����һ����
        tmp = cellfun(@(x)svm_kernel(psi(:,jj), x, kernel_type_ev, cmd_ev), feature{tt},'un',0);
        % K_matΪ n_A*n_feature ����
        K_mat(jj,:) = reshape(cell2mat(tmp), 1, []); 
    end
    % ���aplhaɸѡ��K_mat�Ժ�Ľ�� 
    alphaK_mat  = reshape(alpha'*K_mat, size(flowvar{tt}));
    K = K + sum(sum(alphaK_mat.*flowvar{tt}));
end

    


    
































