% ======================================================================= %
%
% ����� SSVM ѵ���������� 2015.6.17
% ���� active structured learning �е�α�����д��Fig.4��
% ��Ҫ��Ϊ3�����裺
%   1. ���� CXSL_ILP ���㵱ǰw�µ���ѷ��䷽�� z^
%
%   2. �����ݶ� U(x,z*,z^)=fai(x,z*)-fai(x,z^)�����õ�a��b��ֵ����ʽ12��13��
%
%   3. ͨ��a��b�ⷽ��14������º�� w�� kexi������gap�Ĵ�С
%
%   �����з������ڴ治������⣬��Ҫ��û���������ڴ棬���ͨ��Ԥ��������ռ䡢
% ÿ����50�� clear һ�±������������������
% 
%
% ======================================================================= %
clear;close all;
% ���� CXSL_Test_Linear_all �м���õ� w ��Ϊ��ʼֵ��ʵ�ʷ���Ч�������ã�
if 1
    load('E:\datasets\first editon\training datasets\N2DL-HeLa\01_2-16_track\�ṹ��ѧϰ\initial_w_New.mat');
    % ע�� w ��˳������
    w = [ wij, wit, wid, wiv, wmj, wsj ]';
    clear wij wit wid wiv wmj wsj
else
    % ���ѡȡw�����ӿ��� w �ɸ���
    rng(0); 
    w = rand(42,1);
%     w = zeros(42,1);
end
% ������������ N �� ���������е�֡�� frame
N = 1;
frame = 3;
s_frame = zeros(N,1);
e_frame = zeros(N,1);
% Ŀǰ��gt��֡���������ȡ����Ӱ��
gt_frame = 70;

% ѡ��ȡ����ʽ
% 1Ϊ����ȡ����2Ϊ����ȡ����3Ϊ���ȡ��
sample_method = 2;
switch sample_method
    case 1
        % ȡ������1������ȡ��
        % �� 1-5, 5-9, 9-13, 13-17 �����ķ���ȡ����
        for ind=1:N
            s_frame(ind) = (ind - 1)*(frame - 1) + 1;
            e_frame(ind) = s_frame(ind) + frame - 1;
        end

    case 2
        % ȡ������2������ȡ��
        % �� 1-5��2-6��3-7 �����ķ���ȡ����
        for ind=1:N
            s_frame(ind) = ind;
            e_frame(ind) = s_frame(ind) + frame - 1;
        end
        
    case 3
        % ȡ������3�����ȡ��
        rng(0);
        s_frame = randi([1 gt_frame-frame+1], [N 1]);
        e_frame = s_frame + frame - 1;        
end

% ��ӡ������Ϣ
disp(['  ѵ����ѡȡ', num2str(N), '��������']);
for ind=1:N
    disp(['  ', num2str(s_frame(ind)), '����', num2str(e_frame(ind)), '֡...']);
end

% ============================== ȫ�ֱ��� =============================== %
% ȫ���������ն��������������ֲ��������豣��

% ��ʼ���� ����A��B��ѭ���������� iter����϶��ֵ gap
iter = 300;
gap = 0.05; % ���� O(1/gap) �������ٶȣ�Ӧ���ڰ�ѭ���������
gap_cur = zeros(iter,1); % ��¼ÿ�εõ���gap
gap_cur(1) = gap* 2;

time = zeros(iter,1); % ��¼ÿ��ѭ�����õ�ʱ��
sample_loss = zeros(iter,N); % ��¼ÿһ����ÿ����������ʧ����
now_loss = zeros(iter,1); % ��¼ÿһ����������ʧ������ֵ
gamma = zeros(iter,1); % ����gamma

W = cell(iter,1); % W ����ۺ�Ȩֵw
Wi = cell(iter,N); % Wi ���ÿ��ѭ�����ض��������º��Wi
Wi_old = cell(1,N); % Wi_old ���ÿ��ѭ�����ض���������ǰ��Wi
Ws = cell(iter,N);
L = zeros(iter,1); % L ����ۺ�L
Li = zeros(iter,N); % Li ���ÿ��ѡ�е��������º��L
Li_old = zeros(1,N); % Li_old ���ÿ��ѡ�е���������ǰ��L
ls = 1; % ����ƽ����ʧ����

W{1} = w; % ȫ��������W��Ҫ�趨��ֵw
for i=1:N
    Wi{1,i} = w; % ����������WiҲ�趨��ֵΪw������ᵼ�µ�����BCFW��FW֮��Ľ����һ��
    Wi_old{i} = w;
end

% ------------------- ����ΪMP-BCFW���ӵı��� ----------------- %
N_set = 1000; % ��������С����Ϊһ��������
M_pass = 1000; % ���н���oracle����������ʵ��Ϊ��̬���������ﵽ������ǿ��ֹͣ
T_last = 10; % ���T��ѭ���������Ƴ�Wi
Fobj = zeros(iter,N); % ���ÿ��������Ŀ�꺯��

t = 1;
% ======================================================================= %

fai_x_z_hat = cell(N,1);
delta_zstar_zhat = zeros(N,1);

% ============ ����2�ı��� ========== %
U_x_zstar_zhat = cell(N,1);
% ======================================================================== %

disp('  Ԥ����Ŀ�꺯����Լ������...');
%% ��ǰ����� fai(x,z^) fai(x,z*)�͡�(z*,z^)������Լ��������ѭ������װĿ�꺯���������
fij = cell(N,1);
fit = cell(N,1);
fid = cell(N,1);
fiv = cell(N,1);
fmj = cell(N,1);
fsj = cell(N,1);

F = cell(N,1);
for ind=1:N
    F{ind} = lmi;
end

fai_x_z = cell(N,1);
sum_cost = cell(N,1);
fai_x_z_star = cell(N,1);

tic;
% ���� fai(x,z)�͡�(z*,z)����������̱���
for ind=1:N
    % ----------------------------------------- %
    % ������¼����̱�����Ԥ�ȼ���� fai(x,z)�� ��(z*,z)
    [ fij{ind} fit{ind} fid{ind} fiv{ind} fmj{ind} fsj{ind} fai_x_z{ind} sum_cost{ind} ] =...
        CXSL_Calculate_Fai_With_Loss_New( s_frame(ind), e_frame(ind) );
    % ����Լ������ F������ CXSL_Calculate_Constraint_New_Conflict �������
    % �� BundleMethod_Output_Test �е�ͬ������һ��
    % ----------------------------------------- %
    % 2015.7.6 ʹ�����µ�ì��Լ������22ì��Լ����
    [ F{ind} ] = CXSL_Calculate_Constraint_New_Conflict( s_frame(ind), e_frame(ind),...
        fij{ind}, fit{ind}, fid{ind}, fiv{ind}, fmj{ind}, fsj{ind} );
    % ----------------------------------------- %
	% �����׼���е�fai(x,z*)
	[ fai_x_z_star{ind} ] = CXSL_Calculate_fai_x_zstar_New( s_frame(ind), e_frame(ind), 'star'); 
    % ----------------------------------------- %
end
toc;

%% ����ǰѭ������tС�����ޣ���gap������Ҫ��ʱ������ѭ�����㣬�������󾫶Ȼ��������޸�gap��iter�����д�cell����

options = sdpsettings('verbose', 0, 'solver', 'cplex', 'saveduals', 0); % cplex���÷ŵ�ѭ����
rng(0); % �������ѡ�񲿷֣���Ҫ�趨����

while t <= iter && ls*N >= gap

    % ��¼��ÿ��ѭ�����õ�ʱ��
    t_start = clock;
    disp('  ==========================');
    disp(['  ��ʼ�� ', num2str(t), ' ��ѭ��...']);
    
	%% step1. ����һ�ξ�ȷoracle
    tic;
    
    	%% 1. �������w�µľ�ȷoracle
    
    ind = randi(N);
    disp(['      �������� ', num2str(ind), '...']); 
    % ��ʱ�齨Ŀ�꺯�������
    object_function = dot(W{t}, fai_x_z{ind}) + sum_cost{ind};
    sol = solvesdp( F{ind}, -object_function, options );

    % ����õ��ĸ���������ֵ
    if sol.problem == 0      
        fai_x_z_hat{ind} = value(fai_x_z{ind});
        delta_zstar_zhat(ind) = value(sum_cost{ind});
        Fobj(t,ind) = value(object_function); % ��һ�μ���oracle�õ���Ŀ�꺯��ֵ
    else
        sol.info
        yalmiperror(sol.problem)
    end
   
    disp('      ����Ȩ���� w...');
    
        %% 2. ���� ��(x,z*,z^)=fai(x,z*)-fai(x,z^) �ݶ�
    % �� U ������������� ��
    
%     for ind=1:N
        % �ݶ� ������fai(x,z*)-fai(x,z^)
        U_x_zstar_zhat{ind} = fai_x_z_star{ind} - fai_x_z_hat{ind};
        % sum( ��(x,z*,z^) )
        sum_U = U_x_zstar_zhat{ind};
        % ����ÿһ����ÿ����������ʧ����
        sample_loss(t, ind) = delta_zstar_zhat(ind);
%     end
    % =================================================================== %
    % sum( ��(z*,z^) ��
    sum_delta = sum(sample_loss(t,:));
    now_loss(t) = sum_delta; 
    fprintf('      ��ǰ������ʧ������(z*,z^):\t%f\n', sum_delta);

        %% 3. ������Ų���������w
    
    % ����ͷ��� lambda ��
    lambda = 1e-2;
	
    Ws{t,ind} = sum_U/(lambda*N);
    ls = sum_delta/N;
    
    % ����gap��gap��ֵ��������lambda����仯���޷�ȷ����������˻�����loss��gap�ȽϺ���
    gap_cur(t+1) = lambda*(Wi_old{ind}- Ws{t,ind})'*W{t}- Li_old(ind)+ ls;
    % ���㲽��gamma
    gamma(t) = gap_cur(t+1)/(lambda*norm(Wi_old{ind} - Ws{t,ind})^2);
    
    % ���� wi��Li�������º��w������ W{t+1,ind}��
    Wi{t+1,ind} = (1- gamma(t))*Wi_old{ind}+ gamma(t)*Ws{t,ind};
    Li(t+1,ind) = (1- gamma(t))*Li_old(ind)+ gamma(t)*ls;

    % ���� w��L�������º��w������ W{t+1,N+1}��
    W{t+1} = W{t} + Wi{t+1,ind} - Wi_old{ind};
    L(t+1) = L(t) + Li(ind) - Li_old(ind); % ����L�ƺ�ûʲô�ã�
    Wi_old{ind} = Wi{t+1,ind};
    Li_old(ind) = Li(t+1,ind);
    
    fprintf('      ��ż��϶gap:\t%f\n', gap_cur(t+1));

        %% 4. �жϹ����������Ƿ񳬹��趨ֵNset
        N_work = 0;
        for uu=1:t
            if ~isempty(Wi{uu,ind})
                N_work = N_work + 1;
            end
        end
        % ��������������ʱ��Ҫɾȥ���Ծ��һ��Wi
        if N_work>N_set
            %
            % δ�����
            %
        end
      
    %% step2. ������ɴν���oracle
    t_exact = toc;
    
        %% 1. �ڹ��������ҵ���ʹĿ�꺯��ֵ�����Ǹ�W
        tbest = find(Fobj==max(Fobj));
        % ����w
        

        
        
        
    % ������һ��ѭ��
    t = t + 1;  
    % ==================================== %
    % ��¼ʱ��
    t_end  = clock;
    time(t) = etime(t_end, t_start);
    fprintf('      ʱ�仨��:\t%1.2f s\n', time(t)); 
end

% ѭ����ɣ���ӡ��Ϣ
if ls*N <= gap
    disp('  �ҵ��˵�ǰgap�µ����Ž⣬�㷨��ֹ');
    
    % �������ŷ��䷽����w
    t_best = t - 1;
    w_best = W{t-1}; % W{t-1}����ȡ����� gap ֵ���Ǹ�w�������µõ��� W{t} gap����������
    gap_best = ls*N;      
else
    disp('  �ﵽ���ѭ���������㷨��ֹ');
    t = t - 1;
    t_best = find(now_loss==min(now_loss)); %  �ҵ���������ʧ��С���Ǹ�w��Ϊ w_best
    w_best = W{t_best};
    gap_best = now_loss(t_best);
end
% �������w�����ڲ�������֡����
save('E:\datasets\first editon\training datasets\N2DL-HeLa\01_2-16_track\�ṹ��ѧϰ\SSVM_Best_W_New.mat', 'w_best');

figure;plot(now_loss);
% subplot(211);  %aver_loss = aver_loss(1:t); 
% subplot(212); plot(gap_cur); %gap_cur = gap_cur(1:t); 
% fprintf('\tw:\t%f\n', w_best);
fprintf('\n\tt_best:\t%d\n', t_best);
fprintf('\tgap_best:\t%f\n', gap_best);
fprintf('\ttime consumption:\t%0.2f min\n', sum(time)/60);   
% w_for_excel = w_best';









