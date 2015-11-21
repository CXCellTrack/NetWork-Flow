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

% �˷���ÿ������������һ��Wi��ʹ���������µ�ʱ�����Ǵ���һ�ε�Wi=0��ʼ
clear;close all;

[ ~, trackpath ] = getpath( 'training' );
% ���� CXSL_Test_Linear_all �м���õ� w ��Ϊ��ʼֵ��ʵ�ʷ���Ч�������ã�
load([ trackpath, '\�ṹ��ѧϰ\initial_w_New.mat']);
% ע�� w ��˳������
w = [ wij,bij, wit,bit, wid,bid, wiv,biv, wmj,bmj, wsj,bsj ]'; % ԭ��b
if 0
    % ��һ�����ݼ�
%     w = [-24.05027785	0.965959752	-0.209700235	0.023655542	-0.901444678	0.915527485	-0.723055368	0.78127324	-22.34659216	-3.45491283	-1.682414322	-5.355960441	-2.391659001	2.862181421	-7.382944338	8.382838223	1.94377663	-0.451290137	-1.07738777	-4.844423375	-1.122913059	-0.801496889	3.907101647	-11.61160994	3.710115534	0.998335816	4.252699702	0.790594494	1.207125853	3.799458373	1.390618031	5.18991389	1.129864864	0.673380786	-2.076937813	-1.97433464	-1.980221778	-0.051210814	0.597328997	-3.897482158]';
    % �ڶ������ݼ�
%     w = [-18.21315239	2.048551055	0.090611096	-0.07830978	-0.681768441	0.091705287	-0.284558766	0.113666465	-16.77170209	-2.820207584	-1.606735489	-2.170929556	-2.000511632	2.42450433	-4.406444861	10.34098417	1.758814312	-0.819906672	-2.159095585	-4.969572233	1.29607646	0.202318113	3.379177651	-14.25716614	7.109097539	3.366559674	2.10659084	-2.499899814	-3.9849466	2.501397849	1.351917771	4.699879458	-0.525793863	-0.261628668	-4.945586679	-1.846739083	-4.998199696	-0.003648138	1.737582541	-8.35761154]';
    % ���������ݼ�
    w = [-17.69956928	1.61055833	-0.001787913	0.199592206	-0.287005286	-1.058810184	-0.194715268	0.170821175	-15.86943474	-2.459629861	-1.102503407	-1.233255838	-1.209682446	-2.849135318	-5.932152347	3.635597276	3.112043952	0.905772748	-0.112320949	-0.658275898	4.373948249	0.504445114	-0.211398811	-3.571706443	0	0	0	0	0	0	0	0	0	0	-2.873279295	-0.074622812	-1.438747225	-0.381464188	-2.838608426	-5.896985464]';
else
    % initial_w�������
    w = zeros(numel(w), 1);
end

% ������������ N �� ���������е�֡�� frame
N = 5;
frame = 13;
s_frame = zeros(N,1);
e_frame = zeros(N,1);
% Ŀǰ��gt��֡���������ȡ����Ӱ��
gt_frame = 65;

% ѡ��ȡ����ʽ
sample_method = 2;
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
for ind=1:N
    disp(['  ', num2str(s_frame(ind)), '����', num2str(e_frame(ind)), '֡...']);
end

% ȫ���������ն��������������ֲ��������豣��

% ��ʼ���� ����A��B��ѭ���������� iter����϶��ֵ gap
iter = 1000;
gap = 0.0010; % ���� O(1/gap) �������ٶȣ�Ӧ���ڰ�ѭ���������

%% ============================== ȫ�ֱ��� =============================== %
gap_cur = zeros(iter,1); % ��¼ÿ�εõ���gap
gamma = zeros(iter,1); % ����gamma

W = cell(iter,1); % W ����ۺ�Ȩֵw
Wavg = cell(iter,1); % W ����ۺ�Ȩֵw
Wi = cell(iter,N); % Wi�������Ȩֵw
Wavg{1} = w;
W{1} = w; % ȫ��������W��Ҫ�趨��ֵw
for i=1:N
    Wi{1,i} = w; % Wi ���ÿ��ѭ�����ض��������º��Wi
end

L = zeros(iter,1); % L ����ۺ�L
Li = zeros(iter,N);
for i=1:N
    Li(1,i) = 0;
end

t = 0;
time = zeros(iter,1); % ��¼ÿ��ѭ�����õ�ʱ��
sample_loss = zeros(iter,N); % ��¼ÿһ����ÿ����������ʧ����
aver_loss = zeros(iter,1); % ��¼ÿһ����������ʧ������ֵ
% ======================================================================= %
phi_x_z_hat = cell(N,1);
delta_zstar_zhat = zeros(N,1);
% U_x_zstar_zhat = cell(N,1);
% ======================================================================== %
% ѭ����ⲿ�ֲ�������
options = sdpsettings('verbose', 0, 'solver', 'gurobi','cachesolvers',1); % cplex���÷ŵ�ѭ����
rng(0); % �������ѡ�񲿷֣���Ҫ�趨����
random = 1; % ����random��Ϊһ��flag��Ϊ1ʱ�����������Ϊ0ʱ�ǰ�˳�����
ind = 0;

disp('  Ԥ����Ŀ�꺯����Լ������...');

%% ��ǰ����� phi(x,z^) phi(x,z*)�͡�(z*,z^)������Լ��������ѭ������װĿ�꺯���������
use_op_cons = [3 5];

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
sum_cost_all = cell(N,1);
phi_x_z_star = cell(N,1);

tic;
% ���� phi(x,z)�͡�(z*,z)����������̱���
for ii=1:N
    disp('  ==========================');
    disp(['  Ԥ��������',num2str(ii),'��ѵ������...']);
    % ----------------------------------------- %
    % ������¼����̱�����Ԥ�ȼ���� phi(x,z)�� ��(z*,z)
    [ fij{ii} fit{ii} fid{ii} fiv{ii} fmj{ii} fsj{ii} phi_x_z{ii} sum_cost{ii} sum_cost_all{ii} ] =...
        CXSL_Calculate_phi_And_Loss( w, s_frame(ii), e_frame(ii) );
    % ����Լ������ F������ CXSL_Calculate_Constraint_New_Conflict �������
    % �� BundleMethod_Output_Test �е�ͬ������һ��
    % ----------------------------------------- %
    % 2015.7.6 ʹ�����µ�ì��Լ������22ì��Լ����
    [ F{ii} ] = CXSL_Calculate_Constraint_New_Conflict( 'training', use_op_cons, s_frame(ii), e_frame(ii),...
        fij{ii}, fit{ii}, fid{ii}, fiv{ii}, fmj{ii}, fsj{ii} );
    % ----------------------------------------- %
	% �����׼���е�phi(x,z*)
	[ phi_x_z_star{ii} ] = CXSL_Calculate_phi_x_zstar_New( w, s_frame(ii), e_frame(ii), 'star'); 
    % ----------------------------------------- %
end
toc;

%% ����ǰѭ������tС�����ޣ���gap������Ҫ��ʱ������ѭ�����㣬�������󾫶Ȼ��������޸�gap��iter�����д�cell����
% ����ͷ��� lambda ��
n4gap = 2; % ÿN*n4gap��ѭ���������gap��n4gap>1��
lambda = 1e-2;
linesearch = 0;
usecostall = 0;
if usecostall
    disp('��ǰѡ�����ʧ�а������龰��');
    sum_cost = sum_cost_all;
end

while t < iter % && ls*N >= gap
    t = t + 1;
    
    % ��¼��ÿ��ѭ�����õ�ʱ��
    t_start = clock;
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
    
    %% 1. ������� w �£���ǰѡ����������ѷ��䣨���̣�10����
    
    disp(['      �������� ', num2str(ind), '...']); 
    % ��ʱ�齨Ŀ�꺯�������
    if mod(t,n4gap*N)~=0 && t~=1
        object_function = dot(Wavg{t}, phi_x_z{ind}) + sum_cost{ind};
        sol = solvesdp( F{ind}, -object_function, options );

        % ����õ��ĸ���������ֵ
        if sol.problem == 0      
            phi_x_z_hat{ind} = value(phi_x_z{ind});
            delta_zstar_zhat(ind) = value(sum_cost{ind});
        else
            sol.info
            yalmiperror(sol.problem)
        end
    else
        disp('      ����һ��ȫ����������ʧ...');
        for ii=1:N
            disp(['        ����',num2str(ii),'...']);
            object_function = dot(Wavg{t}, phi_x_z{ii}) + sum_cost{ii};
            sol = solvesdp( F{ii}, -object_function, options );
            % ����õ��ĸ���������ֵ
            if sol.problem == 0      
                phi_x_z_hat{ii} = value(phi_x_z{ii});
                delta_zstar_zhat(ii) = value(sum_cost{ii});
            else
                sol.info
                yalmiperror(sol.problem)
            end
        end
    end
        
    disp('      ����Ȩ���� w...');
    
    %% 2. ���� ��(x,z*,z^)=phi(x,z*)-phi(x,z^) �ݶ�
    % �� U ������������� ��
    sum_U = 0;
    if mod(t,n4gap*N)~=0 && t~=1
        U_x_zstar_zhat = phi_x_z_star{ind} - phi_x_z_hat{ind};
        sum_U = U_x_zstar_zhat;
        % ����ÿһ����ÿ����������ʧ����
        sample_loss(t, ind) = delta_zstar_zhat(ind);
    else
        % ��������������ƽ����ʧ
        for ii=1:N
            % �ݶ� ������phi(x,z*)-phi(x,z^)
            U_x_zstar_zhat = phi_x_z_star{ii} - phi_x_z_hat{ii};
            % sum( ��(x,z*,z^) )
            sum_U = sum_U + U_x_zstar_zhat;
            % ����ÿһ����ÿ����������ʧ����
            sample_loss(t, ii) = delta_zstar_zhat(ii);
        end
        
    end
    % =================================================================== %
    % sum( ��(z*,z^) ��
    sum_delta = sum(sample_loss(t,:));
    aver_loss(t) = sum_delta;
    fprintf('      ��ǰ������ʧ������(z*,z^):\t%f\n', sample_loss(t,ind));
    if mod(t,n4gap*N)==0 || t==1
        aver_loss(t) = sum_delta/N;
        fprintf('      ��ǰpass��ƽ����ʧ������(z*,z^):\t%f\n', aver_loss(t));
    end

    %% 3. ������Ų���������Wavg
    Ws = sum_U/(lambda*N);
    ls = sum_delta/N;  
    % ����gap��gap��ֵ��������lambda����仯���޷�ȷ����������˻�����loss��gap�ȽϺ���
    if mod(t,n4gap*N)==0 || t==1
        gap_cur(t) = lambda*(W{t}- Ws)'*W{t}- L(t)+ ls;
    end
    % ���㲽��gamma
    if linesearch
        tmp = lambda*(Wi{t,ind}- Ws)'*W{t}- Li(t,ind)+ ls;
        gamma(t) = tmp/(lambda*norm(Wi{t,ind}- Ws)^2);
        gamma(t) = max([0, min([gamma(t),1])]);
    else
        gamma(t) = 2*N/(2*N + t-1);
    end
    % ���� wi��Li�������º��w������ W{t+1,ind}��
    Wi{t+1,ind} = (1- gamma(t))*Wi{t,ind} + gamma(t)*Ws;
    Li(t+1,ind) = (1- gamma(t))*Li(t,ind) + gamma(t)*ls;
    % ���ڴ���û�ֵ�������������Wi��Li������һ����
    sample_not_used = mysetdiff(1:N,ind);
    for jj=sample_not_used
        % jjΪû���õ����������
        Wi{t+1,jj} = Wi{t,jj};
        Li(t+1,jj) = Li(t,jj); % ֱ�ӽ���Wi������һ��
    end
    
    % ���� w��L�������º��w������ W{t+1,N+1}��
    W{t+1} = W{t} + Wi{t+1,ind} - Wi{t,ind};
    L(t+1) = L(t) + Li(t+1,ind) - Li(t,ind); % ����L�ƺ�ûʲô�ã������ã�2015.11.8��
    
    % ����Wavg
    Wavg{t+1} = (t-1)/(t+1)*Wavg{t} + 2/(t+1)*W{t+1};

    fprintf('      ��ż��϶gap:\t%f\n', gap_cur(t));
    % ==================================== %
    % ��¼ʱ��
    time(t) = etime(clock, t_start);
    fprintf('      ʱ�仨��:\t%1.2f s\n', time(t)); 
    
end

% ѭ����ɣ���ӡ��Ϣ
if ls*N <= gap
    disp('  �ҵ��˵�ǰgap�µ����Ž⣬�㷨��ֹ');
    % �������ŷ��䷽����w
    t_best = t;
    w_best = Wavg{t}; % W{t}����ȡ����� gap ֵ���Ǹ�w
    loss_best = ls*N;      
else
    disp('  �ﵽ���ѭ���������㷨��ֹ');
    
    t_best = find(aver_loss==min(aver_loss(N+1:end))) %  �ҵ���������ʧ��С���Ǹ�w��Ϊ w_best
    t_best = t_best(t_best>N);
    t_best = t_best(end);
    % gapС��ѧϰ��ʧ��һ��С��test���Ҳ��һ��С����������ѧϰ��ʧ�Ƚ�׼��
%     t_best = find(gap_cur==min(gap_cur(gap_cur~=0)));
    w_best = Wavg{t_best};
    loss_best = aver_loss(t_best);
end
% �������w�����ڲ�������֡����
save([ trackpath, '\�ṹ��ѧϰ\SSVM_Best_W_New.mat'], 'w_best');

fprintf('\n\tt_best:\t%d\n', t_best);
fprintf('\tloss_best:\t%f\n', loss_best);
fprintf('\ttime consumption:\t%0.2f min\n', sum(time)/60);   
w_for_excel = w_best';

% ����loss��gap
subplot(211)
    plot(aver_loss, '-*');hold on;
    gapt = 0:n4gap*N:iter;gapt(1) = 1;
    gaploss = aver_loss(gapt);
    plot(gapt,gaploss, 'r-o','linewidth',2);hold off;
subplot(212)
    gapc = gap_cur(gapt);
    plot(gapt,gapc, 'r-o','linewidth',2);
% �Եõ����������߽��б���
if 0
    name = 'loss_5_13_cons35_cost1_init0p_noline_rng';
    lossdir = [ trackpath, '\ѵ�������¼\BCFWavg_paper\new_sample_method\'];
    mkdir(lossdir);
    save([lossdir, name, '.mat'], 'time','w','linesearch','use_op_cons','sample_loss','lambda','w_best','Wavg');
    saveas(1, [lossdir, name, '.fig']);
end







