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
% ���� CXSL_Test_Linear_all �м���õ� w ��Ϊ��ʼֵ
load([ trackpath, '\�ṹ��ѧϰ\initial_w_New.mat']);
% ע�� w ��˳������
w = cell(1,6); % ��wҲ���¼�
if 0
    w{1} = wij'; % ����svm�õ���w
    w{2} = wit';
    w{3} = wid';
    w{4} = wiv';
    w{5} = wmj';
    w{6} = wsj';
else
    w{1} = zeros(numel(wij),1); % ʹ��ȫ0w
    w{2} = zeros(numel(wit),1);
    w{3} = zeros(numel(wid),1);
    w{4} = zeros(numel(wiv),1);
    w{5} = zeros(numel(wmj),1);
    w{6} = zeros(numel(wsj),1);
end


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

%% ------------------ ����w���õ��ı��� ---------------------- %
W = cell(iter,6); % W ����ۺ�Ȩֵw
Wi = cell(iter,N); % Wi�������Ȩֵw
for ev=1:6
    W{1,ev} = w{ev};% ȫ��������W��Ҫ�趨��ֵw
end
for i=1:N
    Wi{1,i} = cell(1,6); % Wi ���ÿ��ѭ�����ض��������º��Wi
    for ev=1:6
        Wi{1,i}{ev} = w{ev};% ȫ��������W��Ҫ�趨��ֵw
    end 
end

L = zeros(iter,6); % L ����ۺ�L
Li = cell(iter,N);
for i=1:N
    Li{1,i} = zeros(1,6);
end

% ����Ͳ�ҪL��
ls = 1; % ����ƽ����ʧ����

%% ------------------ alpha�˺������õ��ı��� ----------------- % 
% ����6���¼��ĺ˺��������Լ���������ѡ�˺���Ϊlinear��poly��rbf��sigmoid��
% 6���¼�����Ϊfij��fit��fid��fiv��fmj��fsj
kernel_type = {'linear','linear','rbf','linear','linear','linear'};
cmd = {'-d 2 -g 1','-d 2 -g 1','-d 2 -g 1','-d 2 -g 1','-d 2 -g 1','-d 2 -g 1'};
% ����ͷ��� lambda �ˣ����ڿ���w����������
lambda = 1e2*ones(1,6);
islinear = strncmp(kernel_type, 'linear', 6); % �߼���������ָʾ��Щ�¼�������

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

alpha_all = cell(iter,1); % ��������alpha�����
for tt=1:iter
    alpha_all{tt} = cell(1,6);
end

phi_y_i = cell(1,N);
for ii=1:N
    phi_y_i{ii} = cell(1,6); % phi_y_iΪd*m���󣬵�������phi�ļ���
end

phi_y_all = cell(iter,1);
for tt=1:iter
    phi_y_all{tt} = cell(1,6); % phi_y_allΪ��������phi�ļ���
end

psi_y_all = cell(iter,1);
for tt=1:iter
    psi_y_all{tt} = cell(1,6); % psi_y_allΪ��������psi�ļ���
end

K_phi = cell(6,1);
phi_x_z_hat = cell(6,1);
n_SV = cell(iter,1);
for tt=1:iter
    n_SV{tt} = zeros(N,6);
end

%% ------------------ ѭ�����õ��ı��� -------------------- %
% ======================================================================= %
% ��ʼ���� ����A��B��ѭ���������� iter����϶��ֵ gap
gap_cur = zeros(iter,1); % ��¼ÿ�εõ���gap
gamma = zeros(iter,1); % ����gamma
t = 0;
time = zeros(iter,1); % ��¼ÿ��ѭ�����õ�ʱ��
sample_loss = zeros(iter,N); % ��¼ÿһ����ÿ����������ʧ����
aver_loss = zeros(iter,1); % ��¼ÿһ����������ʧ������ֵ
% ======================================================================= %
% ѭ����ⲿ�ֲ�������
% options = sdpsettings('verbose', 0, 'solver', 'gurobi');
options = sdpsettings('verbose', 0, 'solver', 'cplex', 'saveduals', 0); % cplex���÷ŵ�ѭ����
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

for ii=1:N
    for ev=1:6
        % ��ʼphi_y����Ϊ��׼�𰸣���Ӧ��psiΪ0��wΪȫ0��
        if ~islinear(ev) % ���ԾͲ���ֵ��
            phi_y_i{ii}{ev} = phi_x_z_star{ii}{ev}; 
        end
    end
end

%% ����ǰѭ������tС�����ޣ���gap������Ҫ��ʱ������ѭ�����㣬�������󾫶Ȼ��������޸�gap��iter�����д�cell����

while (t < iter && ls*N >= gap) || t <= N % �������������������������ÿ�������������õ���
    t = t + 1;
    
    % ��¼��ÿ��ѭ�����õ�ʱ��
    tstart = clock;
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
    disp(['      �������� ', num2str(ind), '...']); 
    % ������������alpha��phi_y
    for ii=1:N
        for ev=1:6
            if ~islinear(ev) % ֻ��Է����Ժ�
                phi_y_all{t}{ev} = [ phi_y_all{t}{ev}, phi_y_i{ii}{ev}];
                % ������phi*���Ƶ���phi_y_i��ͬ�ĳ��ȣ�������õ�psi
                psi_y_i = repmat(phi_x_z_star{ii}{ev},1,size(phi_y_i{ii}{ev},2)) - phi_y_i{ii}{ev};
                psi_y_all{t}{ev} = [ psi_y_all{t}{ev}, psi_y_i];

                % ����alpha
                alpha_all{t}{ev} = [ alpha_all{t}{ev}; alpha_i{ii}{ev}];

                % �������������֧�������������ɹ��鿴��
                n_SV{t}(ii,ev) = numel(alpha_i{ii}{ev});
            end
        end
    end
    
    %% 2. ��ʱ�齨Ŀ�꺯��
    % Ŀ�꺯�����ʽ��
    % H(y�� = sigma[ alpha(i)*K(psi(i), phi(y)) ] + L(y)*lambda*N
    % ��y��������Ϊ
    % H(y�� = sigma[ alpha(i)*y(i)*K(psi(i), feature(i)) ] + L(y)*lambda*N
    % ------------------------------------- %
    K_phi_all = 0; % �����Ժ��¼�Ŀ�꺯����
    W_phi_all = 0; % ���Ժ��¼�Ŀ�꺯����
    for ev=1:6
        if ~islinear(ev) % ֻ��Է����Ժ�
            %% ------------- �����Ժ�Ŀ�꺯�� ----------------- %
            % ʹ�÷����Ժ�ʱ����ڻ��������˺�������
            % ���Ұ���ͬ�¼�����ʹ�ò�ͬ��
            psi = psi_y_all{t}{ev};
            alpha = alpha_all{t}{ev};

            % ʹ��info����¼phi2�������Ϣ������kernel���õ�feature������phi2��
            info.ind = ind;
            info.ev = ev;
            info.s_frame = s_frame(ind); % ѵ��������ʼ֡
            info.e_frame = e_frame(ind); % ѵ����������֡
            % ѡ��ˣ�Ŀǰ��2��feature��������1������Ͳ��ӵ�ԭʼ��
            % �����Ҫ�ж���ǰ����������ĸ�����
            K_phi{ev} = ssvm_kernel_train_my(psi, alpha, info, kernel_type{ev}, cmd{ev},...
                        fij, fit, fid, fiv, fmj, fsj);    

            K_phi_all = K_phi_all + K_phi{ev}/(lambda(ev)*N); % Ŀ�꺯�����ʽ��Ҫ�������¼�������

        else
            %% --------------- ���Ժ�Ŀ�꺯�� ---------------- %
            % ����������˵��ֻ�赼��<w,phi>����
            W_phi_all = W_phi_all + dot(W{t,ev}, phi_x_z{ind}{ev});

        end
    end
    object_function = K_phi_all + W_phi_all + sum_cost{ind};
    sol = solvesdp( F{ind}, -object_function, options );

    %% 3. ����õ��ĸ���������ֵ
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
    
    %% 4. �ж�s��Ӧ��phi֮ǰ�Ƿ���ֹ���������alpha
    % ���õ���phi_x_z_hat��֮ǰ���е�phi���бȽ�
    % ʵ���ϣ����2��y����ͬ�����õ���phi�п�����ͬ�����ͨ���Ƚ�phi��ȷ��y�Ƿ�һ�����Ͻ�
    % ��ȷ�ķ���Ӧ��ͨ���Ƚ�y�����У���ֻҪphi��ͬ��alpha��phi���֮��Ľ������һ��������ޱ�������
    for ev=1:6 % ÿ���¼��ֱ����
        %% --------------- ���Ժ˸��·��� ---------------- %
        if islinear(ev)
            PSI_zstar_zhat = phi_x_z_star{ind}{ev} - phi_x_z_hat{ev};
            Ws = PSI_zstar_zhat/(lambda(ev)*N);
            ls = delta_zstar_zhat/N;  
            % ����gap
            gap_cur = lambda(ev)*(Wi{t,ind}{ev}- Ws)'*W{t,ev}- Li{t,ind}(ev)+ ls;
            if 0
                % line-searchѰ�����gamma
                gamma(t) = gap_cur/(lambda(ev)*norm(Wi{t,ind}{ev}- Ws)^2);
                gamma(gamma>1) = 1; gamma(gamma<0) = 0; % clip to 0-1
            else
                % ��ͨ�������㲽��gamma
                gamma(t) = 2*N/(2*N + t-1);
            end
            % ���� wi��Li�������º��w������ W{t+1,ind}��
            Wi{t+1,ind}{ev} = (1- gamma(t))*Wi{t,ind}{ev} + gamma(t)*Ws;
            Li{t+1,ind}(ev) = (1- gamma(t))*Li{t,ind}(ev) + gamma(t)*ls;
            % ���ڴ���û�ֵ�������������Wi��Li������һ����
            sample_not_used = setdiff(1:N, ind);
            for jj=1:numel(sample_not_used)
                inu = sample_not_used(jj); % û���õ����������
                Wi{t+1,inu}{ev} = Wi{t,inu}{ev}; % ֱ�ӽ���Wi������һ��
                Li{t+1,inu}(ev) = Li{t,inu}(ev);
            end

            % ���� w��L�������º��w������ W{t+1,N+1}��
            W{t+1,ev} = W{t,ev} + Wi{t+1,ind}{ev} - Wi{t,ind}{ev};
        end
        
        %% --------------- �����Ժ˸��·��� ---------------- %
        if ~islinear(ev)
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
        
    end

    %% ͳ����ʧ��ʱ�仨�ѵ�����
    % =================================================================== %
    sample_loss(t, ind) = delta_zstar_zhat;
    aver_loss(t) = delta_zstar_zhat; 
    fprintf('      ��ǰ������ʧ������(z*,z^):\t%f\n', aver_loss(t));
    % ��¼ʱ��
    time(t) = etime(clock, tstart);
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
    t_best
    t_best = t_best(end);
end   

sample_id = mod(t_best,N); % ����tbest�õ���ʱ���������
sample_id(sample_id==0) = N;

% �õ����յ�A��alpha���ڲ�����ʹ����2������
w_best = W(t_best,:);
A_best = cell(1,6); % ����A=psi/(��N)
for i=1:ev
    A_best{i} = psi_y_all{t_best}{i}/(lambda(i)*N);
end
alpha_best = alpha_all{t_best};
n_SV_best = n_SV{t_best};
loss_best = aver_loss(t_best);  

% �������w�����ڲ�������֡����
save([ trackpath, '\�ṹ��ѧϰ\SSVM_Best_W_New.mat'], 'w_best', 'A_best','alpha_best','kernel_type','cmd','islinear');

fprintf('\n\tt_best:\t%d\n', t_best);
fprintf('\tgap_best:\t%f\n', loss_best);
fprintf('\ttime consumption:\t%0.2f min\n', sum(time)/60);   

plot(aver_loss, '-*');
% �Եõ����������߽��б���
if 0
    name = 'loss_5_13_y';
    lossdir = [ trackpath, '\ѵ�������¼\BCFW_my_hunhe\'];
    mkdir(lossdir);
    save([lossdir, name, '.mat'], 'w_best', 'A_best','alpha_best','kernel_type','cmd','islinear');
    saveas(1, [lossdir, name, '.fig']);
end







