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

[ ~, trackpath ] = getpath( 'training' );
load([ trackpath, '\�ṹ��ѧϰ\initial_w_New.mat']);
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
for ind=1:N
    disp(['  ', num2str(s_frame(ind)), '����', num2str(e_frame(ind)), '֡...']);
end

% ��ʼ���� ����A��B��ѭ���������� iter����϶��ֵ gap
iter = 50;
gap = 0.0010; % ���� O(1/gap) �������ٶȣ�Ӧ���ڰ�ѭ���������

%% ============================== ȫ�ֱ��� =============================== %
% ȫ���������ն��������������ֲ��������豣��
gap_cur = zeros(iter,1); % ��¼ÿ�εõ���gap

A = cell(iter,1);
B = cell(iter,1);
W = cell(iter,1); % W ���Ȩֵw
W{1} = w;
ls = 1;
t = 0;
time = zeros(iter,1); % ��¼ÿ��ѭ�����õ�ʱ��

sample_loss = zeros(iter,N); % ��¼ÿһ����ÿ����������ʧ����
aver_loss = zeros(iter,1); % ��¼ÿһ����������ʧ������ֵ
% ======================================================================= %
% ����ѧϰ���м�Ͽ��Ļ����Դ����������� A��B��w �ȼ���ѭ��
% ============================== �ֲ����� ================================ %
% ��ÿһ��ѭ������Щ������ֵ���ᱻ�滻����Ԥ����ռ���Լӿ��ٶ�
% ============ ����1�ı��� ========== %
phi_x_z_hat = cell(N,1);
delta_zstar_zhat = zeros(N,1);
psi_zstar_zhat = cell(N,1);
% ======================================================================= %
disp('  Ԥ����Ŀ�꺯����Լ������...');

%% ��ǰ����� fai(x,z^) fai(x,z*)�͡�(z*,z^)������Լ��������ѭ������װĿ�꺯���������
use_op_cons = [3 5];

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

phi_x_z = cell(N,1);
phi_x_z_star = cell(N,1);
sum_cost = cell(N,1);
sum_cost_all = cell(N,1);

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

%% ����ǰѭ������tС�����ޣ���gap������Ҫ��ʱ������ѭ�����㣬�������󾫶Ȼ��������޸�gap��iter�������cell����
options = sdpsettings('verbose', 0, 'solver', 'gurobi','cachesolvers',1); % cplex���÷ŵ�ѭ����
% ����ͷ��� lambda �ˣ����ԽС����wԽ�������ٶ��ƺ���ӿ죩
lambda = 1e-8;
usecostall = 0;
if usecostall
    sum_cost = sum_cost_all;
end

while t < iter && ls >= gap

    t = t + 1;
    % ��¼��ÿ��ѭ�����õ�ʱ��
    tic;
    disp('  ==========================');
    disp(['  ��ʼ�� ', num2str(t), ' ��ѭ��...']);
    
    %% 1. ������� w �£�ÿ����������ѷ��䣨���̣�10����
    
    for ind=1:N
        disp(['      �������� ', num2str(ind), '...']); 
        % ��ʱ�齨Ŀ�꺯�������
        object_function = dot(W{t}, phi_x_z{ind}) + sum_cost{ind};
        sol = solvesdp( F{ind}, -object_function, options );

        % ����õ��ĸ���������ֵ
        if sol.problem == 0
            % ��Щ���̱�����ֵ���Ǳ����
            if 0
%                 for zhen = s_frame(ind):e_frame(ind)-1
%                     fij{ind}{zhen} = round(value(fij{ind}{zhen})) ;
%                     fid{ind}{zhen} = round(value(fid{ind}{zhen})) ;
%                     fiv{ind}{zhen} = round(value(fiv{ind}{zhen})) ;
%                     fit{ind}{zhen} = round(value(fit{ind}{zhen})) ;
%                 end
%                 for zhen = s_frame(ind)+1:e_frame(ind)
%                     fsj{ind}{zhen} = round(value(fsj{ind}{zhen})) ;
%                     fmj{ind}{zhen} = round(value(fmj{ind}{zhen})) ;
%                 end    
%                 objfun(ind) = value(object_function);
            end
            phi_x_z_hat{ind} = value(phi_x_z{ind});
            delta_zstar_zhat(ind) = value(sum_cost{ind});
        else
            sol.info
            yalmiperror(sol.problem)
        end
    end

    disp('      ����Ȩ���� w...');
    
    %% 2. ���� ��(x,z*,z^)=fai(x,z*)-fai(x,z^) �ݶȣ����a��b��ֵ
    % �� U ������������� ��
    for ind=1:N
        % �ݶ� ������fai(x,z*)-fai(x,z^) ���� bmrm �������õ��� fai(x,z^)-fai(x,z*)
        % �����ķ����ƺ������������ô��뷽��
        psi_zstar_zhat{ind} = phi_x_z_hat{ind} - phi_x_z_star{ind};
    end
    
    % ======================================================================= %
    % ���չ�ʽ���� At = -1/N* sum( ��(x,z*,z^) )
    % ���չ�ʽ���� Bt = -1/N* sum( ��(z*,z^)+ w'* ��(x,z*,z^) ) - w'*At
    
    sum_U = zeros(size(W{t}));
    for ind=1:N
        % sum( ��(x,z*,z^) )
        sum_U = sum_U + psi_zstar_zhat{ind};
        % ����ÿһ����ÿ����������ʧ����
        sample_loss(t, ind) = delta_zstar_zhat(ind);
    end
    % sum( ��(z*,z^) ��
    sum_delta = sum(sample_loss(t,:));
    
    A{t} = 1/N * sum_U;
    B{t} = 1/N *( sum_delta + dot(W{t}, sum_U) ) - dot(W{t}, A{t});
    aver_loss(t) = sum_delta/N;
    ls = aver_loss(t);
    
    fprintf('      ƽ����ʧ������(z*,z^):\t%f\n', aver_loss(t));

    %% 3. ͨ����ⷽ�̣�14��������w������һ�����ι滮���⣩
	% �����º��w������ W{t+1}��
    [ kexi W{t+1} obj] = CXSL_Update_W_For_BMRM( A, B, lambda );

    % ����gap
%     lower_bound = zeros(t,1);
%     for i=1:t
%         lower_bound(i) = A{i}'*W{t+1} + B{i}; % ƽ���ɳڱ������½�
%     end
    % ==================================== %
    % gap �Ķ��壿 ����ط�������
    % J����ÿһ�ֵ� �˦�(w) + �ƣ�one slack���������ι滮���Ǹ�Ŀ�꺯����
    % ֱ������ʧ������Ϊgap_cur gap_cur(t+1)��ʾ����t��ѭ�������ʧ������С
    gap_cur(t) = aver_loss(t);
    
    % ==================================== %
    % ��¼ʱ��
    time(t) = toc;
%     fprintf('      ���Ƽ�϶ ��:\t%f\n', gap_cur);
    fprintf('      ʱ�仨��:\t%1.2f s\n', time(t)); 
   
end

% ѭ����ɣ���ӡ��Ϣ
if gap_cur(t) <= gap
    disp('  �ҵ��˵�ǰgap�µ����Ž⣬�㷨��ֹ');
    
    % �������ŷ��䷽����w
    t_best = t;
    w_best = W{t}; % W{t-1}����ȡ����� gap ֵ���Ǹ�w�������µõ��� W{t} gap����������
    loss_best = gap_cur(t);      
else
    disp('  �ﵽ���ѭ���������㷨��ֹ');
    t_best = find(aver_loss==min(aver_loss(N+1:end))) %  �ҵ���������ʧ��С���Ǹ�w��Ϊ w_best
    t_best = t_best(end);
    w_best = W{t_best};
    loss_best = aver_loss(t_best);
end
% �������w�����ڲ�������֡����
save([ trackpath, '\�ṹ��ѧϰ\SSVM_Best_W_New.mat'], 'w_best');

% fprintf('\tw:\t%f\n', w_best);
fprintf('\tt_best:\t%d\n', t_best);
fprintf('\tgap_cur:\t%f\n', loss_best);
fprintf('\ttime consumption:\t%0.2f min\n', sum(time)/60);   
w_for_excel = w_best';

plot(aver_loss, '-*');
% �Եõ����������߽��б���
if 0
    name = 'loss_5_13_initwp_1e-8';
    lossdir = [ trackpath, '\ѵ�������¼\BMRM\'];
    mkdir(lossdir);
    save([lossdir, name, '.mat'], 'aver_loss','sample_loss','w_best','W');
    saveas(1, [lossdir, name, '.fig']);
end








