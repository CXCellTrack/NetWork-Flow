function KKK = ssvm_kernel_train_paper(n_SV, alpha_i, info, kernel_type, cmd,...
                    fij, fit, fid, fiv, fmj, fsj, s_frame, e_frame,...
                    feature_fij, feature_fit, feature_fid, feature_fiv, feature_fmj, feature_fsj)    

load([ trackpath, '\GT\GT_after_hand_tune\GT_Flow_Variables_New.mat']); % �����׼��
% ��ssvmѵ��ʱ��phi1�ǳ�������������phi2�Ǳ���
% ��Ҫ��phi2���feature*y����ʹ�ú˺�����ע��feature��һ������
y_star_feature = cell(info.N,1); % y*��Ӧ������
y_star = cell(info.N,1); % y*

switch info.ev % �����¼�ѡ��feature�����̱���
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
    disp(['  ��������',num2str(ind), '...'])
    % fe1Ϊy*������fe2Ϊ��������
    fe1 = y_star_feature{ind};
    y1 = y_star{ind};
    fe2 = sample_feature;
    y2 = flowvar;
    
    frame_result = binvar(numel(fe1),numel(fe2));
    for ii=1:numel(fe1)
        for jj=1:numel(fe2) % ii,jj��ʾ��ǰ�ߵ�i֡�ͺ��ߵ�j֡
            
            K_mat = zeros(numel(fe1{ii}),numel(fe2{jj}));
            for kk=1:numel(fe1{ii})
                for mm=1:numel(fe2{jj})
                    K_mat(kk,mm) = svm_kernel(fe1{ii}{kk}, fe2{jj}{mm}, 'linear', ''); %kernel_type, cmd
                end
            end
            % ��y1��y2��ֱ
            tmpy1 = reshape(y1{ii}, 1, []);
            tmpy2 = reshape(y2{jj}, [], 1);
            % ����൱�� ��i��j( yi*yj*k(fei,fej) )
            frame_result(ii,jj) = tmpy1*K_mat*tmpy2;
        end
    end
    
    % ����ͻ�Ҫ*��������֧�������ĸ���
    KKK = KKK + n_SV(ind, info.ev)*sum(frame_result(:));

end
toc
% �����������5��ѭ����������ǰ��֡������֡��ǰ����������������������ǳ���ʱ��
% ��<f*, f>����*alpha����Ϊalpha��һ��Ϊ1













% dims = numel(feature{2}{1});
n_A = size(psi,2);
K = 0; % �൱������ʱ<phi1, phi2>��ֵ

for i=1:numel(sample_feature)
    if i==1 && (info.ev==5 || info.ev==6)
        alphaK_mat = [];
    else
        K_mat = zeros(n_A, numel(sample_feature{i})); % ���psi��feature�ĺ˺����Ľ��
        % psi��ÿһ����feature��ÿһ����������
        for jj=1:n_A
            tmp = cellfun(@(x)svm_kernel(psi(:,jj), x, kernel_type, cmd), sample_feature{i},'un',0);
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
    


    
































