% ======================================================================= %
%
% ����� SSVM ѵ���������� 2015.6.17
% ���� active structured learning �е�α�����д��Fig.4��
% ��Ҫ��Ϊ3�����裺
%   1. ���� CXSL_ILP ���㵱ǰw�µ���ѷ��䷽�� z^
%
%   2. �����ݶ� U(x,z*,z^)=phi(x,z*)-phi(x,z^)�����õ�a��b��ֵ����ʽ12��13��
%
%   3. ͨ��a��b�ⷽ��14������º�� w�� kexi������gap�Ĵ�С
%
%   �����з������ڴ治������⣬��Ҫ��û���������ڴ棬���ͨ��Ԥ��������ռ䡢
% ÿ����50�� clear һ�±������������������
% 
%
% ======================================================================= %
clear;close all;

[ ~, trackpath ] = getpath( 'training' );
load([ trackpath, '\�ṹ��ѧϰ\Feature_Plus_New.mat']); % ��������

% ������������ N �� ���������е�֡�� frame
N = 5;
frame = 13;
s_frame = zeros(N,1);
e_frame = zeros(N,1);
% Ŀǰ��gt��֡���������ȡ����Ӱ��
gt_frame = 65;

% ѡ��ȡ����ʽ
% 1Ϊ����ȡ����2Ϊ����ȡ����3Ϊ���ȡ��
sample_method = 1;
switch sample_method
    case 1 % ȡ������1������ȡ��
        % �� 1-5, 6-10, 11-15, 16-20 �����ķ���ȡ����
        for ind=1:N
            s_frame(ind) = (ind - 1)*frame + 1;
            e_frame(ind) = s_frame(ind) + frame - 1;
        end

    case 2 % ȡ������2������ȡ��
        % �� 1-5��2-6��3-7 �����ķ���ȡ����
        for ind=1:N
            s_frame(ind) = ind;
            e_frame(ind) = s_frame(ind) + frame - 1;
        end
       
    case 3 % ȡ������3�����ȡ��
        rng(0);
        s_frame = randi([1 gt_frame-frame+1], [N 1]);
        e_frame = s_frame + frame - 1;        
end

% ��ӡ������Ϣ
disp(['  ѵ����ѡȡ', num2str(N), '��������']);
for ind=1:N
    disp(['  ', num2str(s_frame(ind)), '����', num2str(e_frame(ind)), '֡...']);
end

%% ����ѭ����������iter������ָ��gap
iter = 50;
gap = 0.0010;

%% ------------------ �˺������õ��ı��� ----------------- % 
% ����6���¼��ĺ˺��������Լ���������ѡ�˺���Ϊlinear��poly��rbf��sigmoid��
% 6���¼�����Ϊfij��fit��fid��fiv��fmj��fsj
kernel_type = {'linear','linear','rbf','linear','linear','linear'};
cmd = {'-d 2 -g 1','-d 2 -g 1','-d 2 -g 1','-d 2 -g 1','-d 2 -g 1','-d 2 -g 1'};

alpha_i = cell(1,N); % N��������ÿ������һ��alpha_i
for ii=1:N
    alpha_i{ii} = cell(6,1); % 6��ϸ���¼�
    % alpha�е�ÿ����Ϊ������
end
for ii=1:N
    for ev=1:6
        alpha_i{ii}{ev} = 1;
    end
end

alpha_all = cell(iter,1);
for tt=1:iter
    alpha_all{tt} = cell(6,1);
end

phi_y_i = cell(1,N);
for ii=1:N
    phi_y_i{ii} = cell(6,1); % phi_y_iΪd*m���󣬵�������phi�ļ���
end

phi_y_all = cell(iter,1);
for tt=1:iter
    phi_y_all{tt} = cell(6,1); % phi_y_allΪ��������phi�ļ���
end

K_phi = cell(6,1);
phi_x_z_hat = cell(6,1);

%% ------------------ ѭ�����õ��ı��� -------------------- %
% ======================================================================= %
% ��ʼ���� ����A��B��ѭ���������� iter����϶��ֵ gap
gap_cur = zeros(iter,1); % ��¼ÿ�εõ���gap
gamma = zeros(iter,1); % ����gamma
ls = 1;
t = 0;
time = zeros(iter,1); % ��¼ÿ��ѭ�����õ�ʱ��
sample_loss = zeros(iter,N); % ��¼ÿһ����ÿ����������ʧ����
aver_loss = zeros(iter,1); % ��¼ÿһ����������ʧ������ֵ
% ======================================================================= %
% ѭ����ⲿ�ֲ�������
options = sdpsettings('verbose', 0, 'solver', 'gurobi');
% options = sdpsettings('verbose', 0, 'solver', 'cplex', 'saveduals', 0); % cplex���÷ŵ�ѭ����
rng(0); % �������ѡ�񲿷֣���Ҫ�趨����
random = 0; % ����random��Ϊһ��flag��Ϊ1ʱ�����������Ϊ0ʱ�ǰ�˳�����
ind = 0;
% ======================================================================= %
disp('  Ԥ����Ŀ�꺯����Լ������...');

%% ��ǰ����� phi(x,z^) phi(x,z*)�͡�(z*,z^)������Լ��������ѭ������װĿ�꺯���������
fij = cell(N,1);
fit = cell(N,1);
fid = cell(N,1);
fiv = cell(N,1);
fmj = cell(N,1);
fsj = cell(N,1);

F = cell(N,1);
for ii=1:N
    F{ii} = lmi;
end

phi_x_z = cell(N,1);
sum_cost = cell(N,1);
phi_x_z_star = cell(N,1);

tic;
% ���� phi(x,z)�͡�(z*,z)����������̱���
for ii=1:N
    disp('  ==========================');
    disp(['  Ԥ��������',num2str(ii),'��ѵ������...']);
    % ----------------------------------------- %
    % ������¼����̱�����Ԥ�ȼ���� phi(x,z)�� ��(z*,z)
    [ fij{ii} fit{ii} fid{ii} fiv{ii} fmj{ii} fsj{ii} phi_x_z{ii} sum_cost{ii} ] =...
        CXSL_Calculate_Event_Fai_And_Loss( s_frame(ii), e_frame(ii) );
    % ����Լ������ F������ CXSL_Calculate_Constraint_New_Conflict �������
    % �� BundleMethod_Output_Test �е�ͬ������һ��
    % ----------------------------------------- %
    % 2015.7.6 ʹ�����µ�ì��Լ������22ì��Լ����
    [ F{ii} ] = CXSL_Calculate_Constraint_New_Conflict( 'training', true, s_frame(ii), e_frame(ii),...
        fij{ii}, fit{ii}, fid{ii}, fiv{ii}, fmj{ii}, fsj{ii} );
    % ----------------------------------------- %
	% �����׼���е�phi(x,z*)
	[ phi_x_z_star{ii} ] = CXSL_Calculate_event_fai_x_zstar( s_frame(ii), e_frame(ii), 'star'); 
    % ----------------------------------------- %
end
toc;

%% ����ǰѭ������tС�����ޣ���gap������Ҫ��ʱ������ѭ�����㣬�������󾫶Ȼ��������޸�gap��iter�����д�cell����
for ii=1:N
    for ev=1:6
        % ��ʼphi_y����Ϊ��׼�𰸣���Ӧ��psiΪ0��wΪȫ0��
        phi_y_i{ii}{ev} = phi_x_z_star{ii}{ev}; 
    end
end

while (t < iter && ls >= gap) || t <= N % �������������������������ÿ�������������õ���
    t = t + 1;
    
    % ��¼��ÿ��ѭ�����õ�ʱ��
    tic;
    disp('  ==========================');
    disp(['  ��ʼ�� ', num2str(t), ' ��ѭ��...']);
    % ����ѡ�������һ���������ǰ�˳����
    if random
        ind = randi(N); % ѡ�е�ind��������Ϊѵ������
    else
        ind = ind + 1;
        if ind==N + 1
            ind = 1;
        end
    end
    
    %% 1. ������� aplha �£���ǰѡ����������ѷ��䣨���̣�10����

    % ����ͷ��� lambda ��
    lambda = 1e-2;
    disp(['      �������� ', num2str(ind), '...']); 
    % ������������alpha��phi_y
    for ii=1:N
        for ev=1:6
            phi_y_all{t}{ev} = [ phi_y_all{t}{ev}, phi_y_i{ii}{ev}];
        end
    end
    for ii=1:N
        for ev=1:6
            alpha_all{t}{ev} = [ alpha_all{t}{ev}; alpha_i{ii}{ev}];
        end
    end

    %% 2. ��ʱ�齨Ŀ�꺯�������
    % Ŀ�꺯�����ʽ��
    % H(y�� = sigma[ alpha(i)*K(psi(i), phi(y)) ] + L(y)*lambda*N
    % ��y��������Ϊ
    % H(y�� = sigma[ alpha(i)*y(i)*K(psi(i), feature(i)) ] + L(y)*lambda*N
    % ------------------------------------- %
    K_phi_all = 0;
    for ev=1:6
        % ʹ�÷����Ժ�ʱ����ڻ��������˺�������
        % ���Ұ���ͬ�¼�����ʹ�ò�ͬ��
        psi = repmat(phi_x_z_star{ind}{ev},1,size(phi_y_all{t}{ev},2)) - phi_y_all{t}{ev};
        alpha = alpha_all{t}{ev};
        
        % ʹ��info����¼phi2�������Ϣ������kernel���õ�feature������phi2��
        info.ind = ind;
        info.ev = ev;
        info.s_frame = s_frame(ind);
        info.e_frame = e_frame(ind);
        % ѡ���
    	K_phi{ev} = ssvm_kernel_train(psi, alpha, info, kernel_type{ev}, cmd{ev},...
            fij, fit, fid, fiv, fmj, fsj,...
    feature_fij_p, feature_fit_p, feature_fid_p, feature_fiv_p, feature_fmj_p, feature_fsj_p);

        K_phi_all = K_phi_all + K_phi{ev}; % Ŀ�꺯�����ʽ��Ҫ�������¼�������
    end
    object_function = K_phi_all + sum_cost{ind}*lambda*N;
    sol = solvesdp( F{ind}, -object_function, options );

    % ����õ��ĸ���������ֵ
    if sol.problem == 0      
        for ev=1:6
            phi_x_z_hat{ev} = value(phi_x_z{ind}{ev});
        end
            delta_zstar_zhat = value(sum_cost{ind});
    else
        sol.info
        yalmiperror(sol.problem)
    end
   
    disp('      ���¶�ż���� ��...');
    
    %% 3. �ж�s��Ӧ��phi֮ǰ�Ƿ���ֹ���������alpha
    % ���õ���phi_x_z_hat��֮ǰ���е�phi���бȽ�
    for ev=1:6 % ÿ���¼��ֱ����
        flag_equal = 0;
        n_phi = size(phi_y_i{ind}{ev},2); % ��ǰ����֧�������ĸ���
        for m=1:n_phi
            if isequal(phi_y_i{ind}{ev}(:,m), phi_x_z_hat{ev}) 
                flag_equal = true;
                break;
            end
        end
        % ------------------------------------------- %
        if flag_equal
            % ������ظ�������һ�ֵ�phi����
            phi_y_i{ind}{ev} = phi_y_i{ind}{ev};
            s = m; % �ҳ�s��λ��
            n_phi_new = n_phi; 
        else
            phi_y_i{ind}{ev} = [ phi_y_i{ind}{ev}, phi_x_z_hat{ev}]; % ���򽫴��ֵõ���phi������һ����
            s = size(phi_y_i{ind}{ev},2); % sΪ�³��ֵ�
            n_phi_new = n_phi + 1; % ���º��֧��������Ŀ
        end

        alpha_vector = zeros(n_phi_new,1);
        alpha_vector(1:n_phi) = alpha_i{ind}{ev};
        s_vector = zeros(n_phi_new,1);
        s_vector(s) = 1; 
        % ����alpha_i
        gamma = 2*N/(2*N + t-1);
        alpha_i{ind}{ev} = (1-gamma)*alpha_vector + s_vector*gamma;
    end

    %% ͳ����ʧ��ʱ�仨�ѵ�����
    % =================================================================== %
    sample_loss(t, ind) = delta_zstar_zhat;
    ls = sum(sample_loss(t,:));
    aver_loss(t) = ls; 
    fprintf('      ��ǰ������ʧ������(z*,z^):\t%f\n', aver_loss(t));
    % ��¼ʱ��
    time(t) = toc;
%     fprintf('      ���Ƽ�϶ ��:\t%f\n', gap_cur);
    fprintf('      ʱ�仨��:\t%1.2f s\n', time(t)); 
    
end

% ѭ����ɣ���ӡ��Ϣ
if ls*N <= gap
    disp('  �ҵ��˵�ǰgap�µ����Ž⣬�㷨��ֹ');
    % �������ŷ��䷽����w
    t_best = t;
else
    disp('  �ﵽ���ѭ���������㷨��ֹ');
    t_best = find(aver_loss==min(aver_loss(aver_loss~=0))); %  �ҵ���������ʧ��С���Ǹ�w��Ϊ w_best   
    t_best = t_best(end);
end   

sample_id = mod(t_best,N); % ����tbest�õ���ʱ���������
sample_id(sample_id==0) = N;

% �õ����յ�A��alpha���ڲ�����ʹ����2������
A_best = phi_y_all{t_best};
alpha_best = alpha_all{t_best};
loss_best = aver_loss(t_best);  

% �������w�����ڲ�������֡����
save([ trackpath, '\�ṹ��ѧϰ\SSVM_Best_W_New.mat'], 'A_best','alpha_best','kernel_type','cmd');

fprintf('\n\tt_best:\t%d\n', t_best);
fprintf('\tgap_best:\t%f\n', loss_best);
fprintf('\ttime consumption:\t%0.2f min\n', sum(time)/60);   

plot(aver_loss, '-*');
% �Եõ����������߽��б���
if 0
    name = 'loss_15_10_y';
    lossdir = [ trackpath, '\ѵ�������¼\BCFW_kernel\'];
    save([lossdir, name, '.mat'], 'sample_loss','A_best','alpha_best','kernel_type','cmd');
    saveas(1, [lossdir, name, '.fig']);
end







