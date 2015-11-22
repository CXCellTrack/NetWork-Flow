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
load([ trackpath, '\�ṹ��ѧϰ\initial_w_New.mat']);
% ע�� w ��˳������
w = cell(1,6); % ��wҲ���¼�
% w{1} = [wij,bij]'; % ����svm�õ���w
% w{2} = [wit,bit]';
% w{3} = [wid,bid]';
% w{4} = [wiv,biv]';
% w{5} = [wmj,bmj]';
% w{6} = [wsj,bsj]';
w{1} = wij'; % ����svm�õ���w
w{2} = wit';
w{3} = wid';
w{4} = wiv';
w{5} = wmj';
w{6} = wsj';

if 0 % ʹ��ȫ0w
    w{1} = zeros(size(w{1}));
    w{2} = zeros(size(w{2}));
    w{3} = zeros(size(w{3}));
    w{4} = zeros(size(w{4}));
    w{5} = zeros(size(w{5}));
    w{6} = zeros(size(w{6}));
end

% ������������ N �� ���������е�֡�� frame
N = 5;
frame = 13;
s_frame = zeros(N,1);
e_frame = zeros(N,1);
% Ŀǰ��gt��֡���������ȡ����Ӱ��
gt_frame = 65;

% ѡ��ȡ����ʽ
sample_method = 1;
switch sample_method
    case 1 % ȡ������1������ȡ����ԭ������
        % �� 1-5, 6-10, 12-15, 16-20 �����ķ���ȡ����
        for ii=1:N
            s_frame(ii) = (ii - 1)*frame + 1;
            e_frame(ii) = s_frame(ii) + frame - 1;
        end

    case 2 % ȡ������2���ص�����ȡ����new sample method��
        % �� 1-5, 5-10, 10-15, 15-20 �����ķ���ȡ����
        for ii=1:N
            s_frame(ii) = (ii - 1)*frame;
            e_frame(ii) = s_frame(ii) + frame;
        end
        s_frame(1) = 1;
        
    case 3 % ȡ������3������ȡ��
        % �� 1-5��2-6��3-7 �����ķ���ȡ����
        for ind=1:N
            s_frame(ind) = ind;
            e_frame(ind) = s_frame(ind) + frame - 1;
        end
       
    case 4 % ȡ������4�����ȡ��
        rng(0);
        s_frame = randi([1 gt_frame-frame+1], [N 1]);
        e_frame = s_frame + frame - 1;        
end

% ��ӡ������Ϣ
disp(['  ѵ����ѡȡ', num2str(N), '��������']);
for ii=1:N
    disp(['  ', num2str(s_frame(ii)), '����', num2str(e_frame(ii)), '֡...']);
end

% ����ѭ����������iter������ָ��gap
iter = 50;
gap = 0.0010;

%% ------------------ ����w���õ��ı��� ---------------------- %
W = cell(iter,6); % W ����ۺ�Ȩֵw
Wavg = cell(iter,6);
Wi = cell(iter,N); % Wi�������Ȩֵw
for ev=1:6
    W{1,ev} = w{ev};% ȫ��������W��Ҫ�趨��ֵw
    Wavg{1,ev} = w{ev};
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

%% ------------------ �˺������õ��ı��� ----------------- % 
% ����6���¼��ĺ˺��������Լ���������ѡ�˺���Ϊlinear��poly��rbf��sigmoid��
% 6���¼�����Ϊfij��fit��fid��fiv��fmj��fsj
kernel_type = {'linear','linear','rbf','linear','linear','linear'};
cmd = {'','','-g 0.1','','',''};
% ����ͷ��� lambda �ˣ����ڿ���w����������
lambda = 1e-2*ones(1,6);
islinear = strncmp(kernel_type, 'linear', 6); % �߼���������ָʾ��Щ�¼�������

alpha_i = cell(iter,1); % N��������ÿ������һ��alpha_i
alpha_avg = cell(iter,1);
for tt=1:iter
    alpha_i{tt} = cell(N,6);
    alpha_avg{tt} = cell(N,6);
    for ii=1:N*6
        alpha_i{tt}(ii) = {1};
        alpha_avg{tt}(ii) = {1}; 
    end
end

phi_y_i = cell(N,6);

y_i = cell(iter,1); % ֧������phi��Ӧ��y
for tt=1:iter
    y_i{tt} = cell(N,6);
end


Y_i_star_K_Ytrain = cell(1,6);
Y_i_K_Ytrain = cell(1,6);
phi_x_z_hat = cell(1,6);
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
options = sdpsettings('verbose', 0, 'solver', 'gurobi');
rng(0); % �������ѡ�񲿷֣���Ҫ�趨����
random = 0; % ����random��Ϊһ��flag��Ϊ1ʱ�����������Ϊ0ʱ�ǰ�˳�����
ind = 0;
% ======================================================================= %
disp('  Ԥ����Ŀ�꺯����Լ������...');

%% ��ǰ����� phi(x,z^) phi(x,z*)�͡�(z*,z^)������Լ��������ѭ������װĿ�꺯���������
use_op_cons = [3 5];

fij = cell(N,1);
fit = cell(N,1);
fid = cell(N,1);
fiv = cell(N,1);
fmj = cell(N,1);
fsj = cell(N,1);
% ���ѭ�������õ���y
Fij = cell(e_frame(N)-1,1);
Fid = cell(e_frame(N)-1,1);
Fiv = cell(e_frame(N)-1,1);
Fit = cell(e_frame(N)-1,1);
Fsj = cell(e_frame(N),1);
Fmj = cell(e_frame(N),1);

F = cell(N,1);
Fkernel = cell(N,1);
for ii=1:N
    F{ii} = lmi;
    Fkernel{ii} = lmi;
end

phi_x_z = cell(N,1);
sum_cost = cell(N,1);
sum_cost_all = cell(N,1);
phi_x_z_star = cell(N,1);

% ���� phi(x,z)�͡�(z*,z)����������̱���
for ii=1:N
    disp('  ==========================');
    disp(['  Ԥ��������',num2str(ii),'��ѵ������...']);
    % ----------------------------------------- %
    % ������¼����̱�����Ԥ�ȼ���� phi(x,z)�� ��(z*,z)
    [ fij{ii} fit{ii} fid{ii} fiv{ii} fmj{ii} fsj{ii} phi_x_z{ii} sum_cost{ii} sum_cost_all{ii} ] =...
        CXSL_Calculate_Event_Fai_And_Loss( s_frame(ii), e_frame(ii) );
    % ----------------------------------------- %
    % 2015.7.6 ʹ�����µ�ì��Լ������22ì��Լ����
    [ F{ii} ] = CXSL_Calculate_Constraint_New_Conflict( 'training', use_op_cons, s_frame(ii), e_frame(ii),...
        fij{ii}, fit{ii}, fid{ii}, fiv{ii}, fmj{ii}, fsj{ii} );
    % ��������е�Լ����ֻ��Ϊ�˼�С������Ĺ�ģ��
    [ Fkernel{ii} ] = calculate_kernel_constraint( islinear, s_frame(ii), e_frame(ii),...
        fij{ii}, fit{ii}, fid{ii}, fiv{ii}, fmj{ii}, fsj{ii} );
    % 2Լ���ϲ�
    if 1
        disp('  ������kernelԼ����ʹ�𰸵�ȫ0��ǿ��Ϊ0');
        F{ii} = [ F{ii}, Fkernel{ii} ];
    end
    % ----------------------------------------- %
	% �����׼���е�phi(x,z*)
	[ phi_x_z_star{ii} ] = CXSL_Calculate_event_fai_x_zstar( s_frame(ii), e_frame(ii), 'star'); 
    % ----------------------------------------- %
end

%% ���ʱ�Ĳ��������������������������ĺ�
use_distrabute = 0;
if ~use_distrabute
    kernel_path = [trackpath, '\��ѵ��\kernel_ff_all.mat'];
    if ~exist(kernel_path, 'file')   
        kernel_ff_all = cell(1,6);
        for ev=1:6
            kernel_ff_all{ev} = cell(gt_frame-1);
        end
        s1 = 1; % fe1�Ŀ�ʼ֡
        e1 = 65; % fe1�Ľ���֡
        s2 = 1; % fe2�Ŀ�ʼ֡
        e2 = 65; % fe2�Ľ���֡
        for ev=1:6 
            if ~islinear(ev)
                kernel_ff_all{ev} = ssvm_pre_cal_all_kernel_paper(kernel_ff_all{ev}, gt_frame, ev, kernel_type{ev}, cmd{ev}, [s1 e1 s2 e2]);
                disp('����kernel����...');tic % ������ô�������Ҳ�ܻ�ʱ��
                save(kernel_path, 'kernel_ff_all','-v7.3');toc
            end
        end
    else
        disp('�������ȼ���õĺ˺���...');
        tic
        load(kernel_path); % ������ֱ�����뼴��
        toc
    end

else
    % ʹ�÷�ɢ�洢
    kernel_ff_all = cell(1,6);
    s1 = 1; % fe1�Ŀ�ʼ֡
    e1 = 80; % fe1�Ľ���֡
    s2 = 1; % fe2�Ŀ�ʼ֡
    e2 = 80; % fe2�Ľ���֡
    for ev=1:6 
        ev_kernel_name = [trackpath, '\��ѵ��\dis_kernel_data\row_1_ev',num2str(ev),'.mat'];
        if ~islinear(ev) && ~exist(ev_kernel_name, 'file') 
            ssvm_pre_cal_all_kernel_paper_dis(gt_frame, ev, kernel_type{ev}, cmd{ev}, [s1 e1 s2 e2]);
            for i_file=s1:e1-1
                % ����ֲ�ʽ�洢mat������ϳɴ������������
                tmp_path = [trackpath, '\��ѵ��\dis_kernel_data\row_',num2str(i_file),'_ev',num2str(ev),'.mat'];
                disp('�����',num2str(i_file),'��kernel����...');
                load(tmp_path);
                kernel_ff_all{ev} = [kernel_ff_all{ev}; rowkernel];
            end
        end
    end
end
    
%% ��y_i��phi_y_i�����ֵ
for ii=1:N
    for ev=1:6
        if ~islinear(ev) % ֻ��Է����Ժ�
            % ע��phi_y_i��i��A������������������������
            % ��ʼphi_y����Ϊ��׼�𰸣���Ӧ��psiΪ0��wΪȫ0��
            phi_y_i{ii,ev} = phi_x_z_star{ii}{ev}; 
            % ʹ�ú�������y_i�����ֵ
            y_i{1}{ii,ev}{1} = init_assgin_y_i( trackpath, s_frame(ii), e_frame(ii), ev);
        end
    end
end

% global Kernel_ev; % ����ȫ�ֱ����Լӿ��ٶ�(�۲��ٶ����������������ǲ�����)

%% ����ǰѭ������tС�����ޣ���gap������Ҫ��ʱ������ѭ�����㣬�������󾫶Ȼ��������޸�gap��iter�����д�cell����
usecostall = 0;
linesearch = 1;
if usecostall
    disp('��ǰѡ�����ʧ�а������龰��');
    sample_cost = sum_cost_all;
else
    sample_cost = sum_cost;
end

while t < iter %&& ls*N >= gap) || t <= N % �������������������������ÿ�������������õ���
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
            % �������������֧�������������ɹ��鿴��
            n_SV{t}(ii,ev) = numel(alpha_i{t}{ii,ev});
        end
    end

    %% 2. ��ʱ�齨Ŀ�꺯�������
    % Ŀ�꺯�����ʽ��
    K_OBJ = 0;
    W_phi_all = 0;
    tic
    for ev=1:6
        if ~islinear(ev) % ֻ��Է����Ժ�
            %% ------------- �����Ժ�Ŀ�꺯�� ----------------- %
            % ʹ�÷����Ժ�ʱ����ڻ��������˺�������
            % ���Ұ���ͬ�¼�����ʹ�ò�ͬ��
%             profile on
            Kernel_ev = kernel_ff_all{ev};
            y_ev = y_i{t}(:,ev); % ȡ�����¼���y_i
            alpha_ev = alpha_avg{t}(:,ev); % ȡ�����¼���alpha_i
            ind_train = ind;
            % 1������<phi_i*,phi>
            % ����Y_i'*K*Y 
            % 2������<phi_i,phi>
            % ����Y_i'*K*Y
            [ Y_i_star_K_Ytrain{ev} Y_i_K_Ytrain{ev} ] = ssvm_kernel_train_paper(y_ev, alpha_ev, Kernel_ev, ev, N, ind_train, s_frame, e_frame,...
                                            fij, fit, fid, fiv, fmj, fsj); 
            % 3���ϳ�Ŀ�꺯��
            K_OBJ = K_OBJ + (Y_i_star_K_Ytrain{ev} - Y_i_K_Ytrain{ev})/(lambda(ev)*N); % Ŀ�꺯�����ʽ��Ҫ�������¼�������
%             profile viewer
            
        else
            %% --------------- ���Ժ�Ŀ�꺯�� ---------------- %
            % ����������˵��ֻ�赼��<w,phi>����
            W_phi_all = W_phi_all + dot(Wavg{t,ev}, phi_x_z{ind}{ev});
        
        end
    end
    toc;
    object_function = K_OBJ + W_phi_all + sample_cost{ind};
    disp('    Ŀ�꺯���齨���');
    sol = solvesdp( F{ind}, -object_function, options );

    %% 3. ����õ��ĸ���������ֵ����y_i������Ϊ��֧��������
    if sol.problem == 0      
        for ev=1:6
            phi_x_z_hat{ev} = value(phi_x_z{ind}{ev});
        end
        delta_zstar_zhat = value(sample_cost{ind});
        
        % ��Ҫ��phi_x_z_hat��Ӧ���Ǹ�YҲ����
        for ff = s_frame(ind):e_frame(ind)-1
            Fij{ff} = round(value(fij{ind}{ff})) ;
            Fid{ff} = round(value(fid{ind}{ff})) ;
            Fiv{ff} = round(value(fiv{ind}{ff})) ;
            Fit{ff} = round(value(fit{ind}{ff})) ;
        end
        for ff = s_frame(ind)+1:e_frame(ind)
            Fsj{ff} = round(value(fsj{ind}{ff})) ;
            Fmj{ff} = round(value(fmj{ind}{ff})) ;
        end
        Fsj(s_frame(ind)) = []; % ȥ����һλ�����ֺ������¼�һ�����ȣ���ʵ������ǰ����һ֡��
        Fmj(s_frame(ind)) = []; % ȥ����һλ�����ֺ������¼�һ�����ȣ���ʵ������ǰ����һ֡��
        
    else
        sol.info
        yalmiperror(sol.problem)
    end
   
    disp('      ���¶�ż���� ��...');
    
    %% 4. �ж�s��Ӧ��phi֮ǰ�Ƿ���ֹ���������alpha_i��y_i
    % ���õ���phi_x_z_hat��֮ǰ���е�phi���бȽ�
    % ʵ���ϣ����2��y����ͬ�����õ���phi�п�����ͬ�����ͨ���Ƚ�phi��ȷ��y�Ƿ�һ�����Ͻ�
    % ��ȷ�ķ���Ӧ��ͨ���Ƚ�y�����У���ֻҪphi��ͬ��alpha��phi���֮��Ľ������һ��������ޱ�������
    %
    % ����������������һ�ִ��������ٸ���ind����
    alpha_i{t+1} = alpha_i{t};
    alpha_avg{t+1} = alpha_avg{t}; % ��2�����д���¼�ѭ����������
    %
    
    for ev=1:6 % ÿ���¼��ֱ����
        %% --------------- ���Ժ˸��·��� ---------------- %
        if islinear(ev)
            PSI_zstar_zhat = phi_x_z_star{ind}{ev} - phi_x_z_hat{ev};
            Ws = PSI_zstar_zhat/(lambda(ev)*N);
            ls = delta_zstar_zhat/N;  
            
            if linesearch
                % line-searchѰ�����gamma
                tmp = lambda(ev)*(Wi{t,ind}{ev}- Ws)'*W{t,ev}- Li{t,ind}(ev)+ ls;
                gamma(t) = tmp/(lambda(ev)*norm(Wi{t,ind}{ev}- Ws)^2);
                gamma(t) = max([0, min([gamma(t),1])]); % clip to 0-1
            else
                % ��ͨ�������㲽��gamma
                gamma(t) = 2*N/(2*N + t-1);
            end
            % ���� wi��Li�������º��w������ W{t+1,ind}��
            Wi{t+1,ind}{ev} = (1- gamma(t))*Wi{t,ind}{ev} + gamma(t)*Ws;
            Li{t+1,ind}(ev) = (1- gamma(t))*Li{t,ind}(ev) + gamma(t)*ls;
            % ���ڴ���û�ֵ�������������Wi��Li������һ����
            sample_not_used = mysetdiff(1:N, ind);
            for jj=sample_not_used
                % jj û���õ����������
                Wi{t+1,jj}{ev} = Wi{t,jj}{ev}; % ֱ�ӽ���Wi������һ��
                Li{t+1,jj}(ev) = Li{t,jj}(ev);
            end

            % ���� w��L�������º��w������ W{t+1,N+1}��
            W{t+1,ev} = W{t,ev} + Wi{t+1,ind}{ev} - Wi{t,ind}{ev};
            
            % ����Wavg
            Wavg{t+1,ev} = (t-1)/(t+1)*Wavg{t,ev} + 2/(t+1)*W{t+1,ev};
        end
        
        %% --------------- �����Ժ˸��·��� ---------------- %
        if ~islinear(ev)
            flag_equal = 0;
            n_phi = size(phi_y_i{ind,ev},2); % ��ǰ����֧�������ĸ���
            for m=1:n_phi
                if isequal(phi_y_i{ind,ev}(:,m), phi_x_z_hat{ev}) 
                    flag_equal = true;
                    break;
                end
            end
            % ------------------------------------------- %
            if flag_equal
                % ������ظ�������һ�ֵ�phi����
                phi_y_i{ind,ev} = phi_y_i{ind,ev};
                s = m; % �ҳ�s��λ��
                n_phi_new = n_phi; 
                y_i{t+1} = y_i{t}; % y_i����
            else
                phi_y_i{ind,ev} = [ phi_y_i{ind,ev}, phi_x_z_hat{ev}]; % ���򽫴��ֵõ���phi������һ����
                s = size(phi_y_i{ind,ev},2); % sΪ�³��ֵ�
                n_phi_new = n_phi + 1; % ���º��֧��������Ŀ
                % ͬʱ���Ӧ��y^ҲҪ���뵽y_i��
                y_i{t+1} = y_i{t}; % ����������������һ�ִ��������ٸ���ind����
                y_i{t+1}{ind,ev}{n_phi_new} = assign_y_hat_into_y_i(ev, ind, s_frame, e_frame,Fij,Fit,Fid,Fiv,Fmj,Fsj);
                
            end

            alpha_vector = zeros(n_phi_new,1);
            alpha_vector(1:n_phi) = alpha_i{t}{ind,ev}; % ��aplha����
            s_vector = zeros(n_phi_new,1);
            s_vector(s) = 1; % ��aplha����
            % ����alpha_i
            gamma_alpha = 2*N/(2*N + t-1);
            alpha_i{t+1}{ind,ev} = (1-gamma_alpha)*alpha_vector + s_vector*gamma_alpha;

            % ����alpha_avg���������Ȳ���������һλ����Ϊ�ȳ���
            if numel(alpha_avg{t}{ind,ev})~=numel(alpha_i{t+1}{ind,ev})
                alpha_avg{t}{ind,ev} = [alpha_avg{t}{ind,ev}; 0];
            end
            alpha_avg{t+1}{ind,ev} = (t-1)/(t+1)*alpha_avg{t}{ind,ev} + 2/(t+1)*alpha_i{t+1}{ind,ev};
            
        end
        
    end

    %% 5. ͳ����ʧ��ʱ�仨�ѵ�����
    % =================================================================== %
    sample_loss(t, ind) = delta_zstar_zhat;
    aver_loss(t) = delta_zstar_zhat; 
    fprintf('      ��ǰ������ʧ������(z*,z^):\t%f\n', aver_loss(t));
    % ��¼ʱ��
    time(t) = etime(clock, tstart);
    fprintf('      ʱ�仨��:\t%1.2f s\n', time(t)); 
    
end

% ѭ����ɣ���ӡ��Ϣ
if ls*N <= gap
    disp('  �ҵ��˵�ǰgap�µ����Ž⣬�㷨��ֹ');
    % �������ŷ��䷽����w
    t_best = t;
else
    disp('  �ﵽ���ѭ���������㷨��ֹ');
    t_best = find(aver_loss==min(aver_loss(N+1:end))) %  �ҵ���������ʧ��С���Ǹ�w��Ϊ w_best
    t_best = t_best(t_best>N);
    t_best = t_best(end);
end   

sample_id = mod(t_best,N); % ����tbest�õ���ʱ���������
sample_id(sample_id==0) = N;

%% �õ����յ�y^��y*��alpha_i���ڲ�����ʹ����2������
w_best = Wavg(t_best,:);
yhat_best = y_i{t_best};

load([ trackpath, '\�ṹ��ѧϰ\Feature_Plus_New.mat']); % ��������
load([ trackpath, '\GT\GT_after_hand_tune\GT_Flow_Variables_New.mat']); % �����׼��
ystar_best = cell(N,6);
feature_best = cell(N,6);
for ii=1:N % ���׼��y*
    ystar_best{ii,1} = Fij(s_frame(ii):e_frame(ii)-1);
    ystar_best{ii,2} = Fit(s_frame(ii):e_frame(ii)-1);
    ystar_best{ii,3} = Fid(s_frame(ii):e_frame(ii)-1);
    ystar_best{ii,4} = Fiv(s_frame(ii):e_frame(ii)-1);
    ystar_best{ii,5} = Fmj(s_frame(ii)+1:e_frame(ii));
    ystar_best{ii,6} = Fsj(s_frame(ii)+1:e_frame(ii));
    
    % �ҳ���Ӧ������
    feature_best{ii,1} = feature_fij_p(s_frame(ii):e_frame(ii)-1);
    feature_best{ii,2} = feature_fit_p(s_frame(ii):e_frame(ii)-1);
    feature_best{ii,3} = feature_fid_p(s_frame(ii):e_frame(ii)-1);
    feature_best{ii,4} = feature_fiv_p(s_frame(ii):e_frame(ii)-1);
    feature_best{ii,5} = feature_fmj_p(s_frame(ii)+1:e_frame(ii));
    feature_best{ii,6} = feature_fsj_p(s_frame(ii)+1:e_frame(ii));
end
% ����һ�����ݵĹ�ģ
for ev=1:6 
    if islinear(ev)
        ystar_best(:,ev) = {[]};
        feature_best(:,ev) = {[]};
    end
end

alpha_best = alpha_avg{t_best};
n_SV_best = n_SV{t_best};

% ����ʹ�� y*-y^ ����Ԥ��
delta_y_best = cell(N,6);
for nn=1:N*6
    delta_y_best{nn} = {};
end
for ev=1:6 % �¼�ѭ����ֻ�к�
    if islinear(ev)
        continue;
    end
    for ii=1:N % ����ѭ��
        for nsv=1:n_SV_best(ii,ev) % ÿ�������е�֧��������Ŀ
            for tt=1:numel(ystar_best{ii,ev}) % ÿ������֡��
                % ��y*-y^������ʱ��2�����Ƕ�ֵ�������ȼ���
                delta_y_best{ii,ev}{nsv,tt} = ystar_best{ii,ev}{tt} - yhat_best{ii,ev}{nsv}{tt};
                % ��ii��������ev���¼��У���nsv��alpha(ii)��Ӧ�ĵ�tt֡
            end
        end
    end
end

loss_best = aver_loss(t_best);  

%% �������w�����ڲ�������֡����

fprintf('\n\tt_best:\t%d\n', t_best);
fprintf('\tgap_best:\t%f\n', loss_best);
fprintf('\ttime consumption:\t%0.2f min\n', sum(time)/60);   

plot(aver_loss, '-*');
% �Եõ����������߽��б���
if 0
    name = 'loss_5_13_initwp_line';
    lossdir = [ trackpath, '\ѵ�������¼\�˼�¼\'];
    mkdir(lossdir);
    save([lossdir, name, '.mat'], 'time','sample_loss','w_best','linesearch',...
        'delta_y_best','feature_best','alpha_best','kernel_type','cmd','lambda','islinear');
    saveas(1, [lossdir, name, '.fig']);
end







