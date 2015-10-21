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
% ���� CXSL_Test_Linear_all �м���õ� w ��Ϊ��ʼֵ��ʵ�ʷ���Ч�������ã�
if 1
    load([ trackpath, '\�ṹ��ѧϰ\initial_w_New.mat']);
    % ע�� w ��˳������
    w = [ wij, wit, wid, wiv, wmj, wsj ]';
    clear wij wit wid wiv wmj wsj;
else
    % ���ѡȡw�����ӿ��� w �ɸ���
    rng(0); 
%     w = rand(40,1);
    w = zeros(40,1);
end
% ������������ N �� ���������е�֡�� frame
N = 1;
frame = 5;
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

% ============================== ȫ�ֱ��� =============================== %
% ȫ���������ն��������������ֲ��������豣��

% ��ʼ���� ����A��B��ѭ���������� iter����϶��ֵ gap
iter = 50;
gap = 0.0010; % ���� O(1/gap) �������ٶȣ�Ӧ���ڰ�ѭ���������
gap_cur = zeros(iter,1); % ��¼ÿ�εõ���gap
gamma = zeros(iter,1); % ����gamma

W = cell(iter,1); % W ����ۺ�Ȩֵw
Wi = w; % Wi ���ÿ��ѭ�����ض��������º��Wi
Wi_old = 0; % Wi_old ���ÿ��ѭ�����ض���������ǰ��Wi
Wavg = cell(iter,1); % Wavg
Ws = 0;

L = zeros(iter,1); % L ����ۺ�L
Li = 0; % Li ���ÿ��ѡ�е��������º��L
Li_old = 0; % Li_old ���ÿ��ѡ�е���������ǰ��L
ls = 1; % ����ƽ����ʧ����

W{1} = w; % ȫ��������W��Ҫ�趨��ֵw
Wavg{1} = w; % Wavg���ó�ֵҲΪw

t = 0;
time = zeros(iter,1); % ��¼ÿ��ѭ�����õ�ʱ��
sample_loss = zeros(iter,N); % ��¼ÿһ����ÿ����������ʧ����
aver_loss = zeros(iter,1); % ��¼ÿһ����������ʧ������ֵ
% ======================================================================= %
% ѭ����ⲿ�ֲ�������
options = sdpsettings('verbose', 0, 'solver', 'cplex', 'saveduals', 0); % cplex���÷ŵ�ѭ����
rng(0); % �������ѡ�񲿷֣���Ҫ�趨����
random = 0; % ����random��Ϊһ��flag��Ϊ1ʱ�����������Ϊ0ʱ�ǰ�˳�����
ind = 0;

disp('  Ԥ����Ŀ�꺯����Լ������...');
%% ��ǰ����� fai(x,z^) fai(x,z*)�͡�(z*,z^)������Լ��������ѭ������װĿ�꺯���������
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

fai_x_z = cell(N,1);
sum_cost = cell(N,1);
fai_x_z_star = cell(N,1);

tic;
% ���� fai(x,z)�͡�(z*,z)����������̱���
for ii=1:N
    disp('  ==========================');
    disp(['  Ԥ��������',num2str(ii),'��ѵ������...']);
    % ----------------------------------------- %
    % ������¼����̱�����Ԥ�ȼ���� fai(x,z)�� ��(z*,z)
    [ fij{ii} fit{ii} fid{ii} fiv{ii} fmj{ii} fsj{ii} fai_x_z{ii} sum_cost{ii} ] =...
        CXSL_Calculate_Fai_And_Loss( s_frame(ii), e_frame(ii) );
    % ����Լ������ F������ CXSL_Calculate_Constraint_New_Conflict �������
    % �� BundleMethod_Output_Test �е�ͬ������һ��
    % ----------------------------------------- %
    % 2015.7.6 ʹ�����µ�ì��Լ������22ì��Լ����
    [ F{ii} ] = CXSL_Calculate_Constraint_New_Conflict( 'training', true, s_frame(ii), e_frame(ii),...
        fij{ii}, fit{ii}, fid{ii}, fiv{ii}, fmj{ii}, fsj{ii} );
    % ----------------------------------------- %
	% �����׼���е�fai(x,z*)
	[ fai_x_z_star{ii} ] = CXSL_Calculate_fai_x_zstar_New( s_frame(ii), e_frame(ii), 'star'); 
    % ----------------------------------------- %
end
toc;

%% ����ǰѭ������tС�����ޣ���gap������Ҫ��ʱ������ѭ�����㣬�������󾫶Ȼ��������޸�gap��iter�����д�cell����

while (t < iter && ls*N >= gap) || t <= N % �������������������������ÿ�������������õ���
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
    
    %% 1. ������� w �£���ǰѡ����������ѷ��䣨���̣�10����
    
    disp(['      �������� ', num2str(ind), '...']); 
    % ��ʱ�齨Ŀ�꺯�������
    object_function = dot(Wavg{t}, fai_x_z{ind}) + sum_cost{ind};
    sol = solvesdp( F{ind}, -object_function, options );

    % ����õ��ĸ���������ֵ
    if sol.problem == 0      
        fai_x_z_hat = value(fai_x_z{ind});
        delta_zstar_zhat = value(sum_cost{ind});
    else
        sol.info
        yalmiperror(sol.problem)
    end
   
    disp('      ����Ȩ���� w...');
    
    %% 2. ���� ��(x,z*,z^)=fai(x,z*)-fai(x,z^) �ݶ�
    % �� U ������������� ��
    
    % �ݶ� ������fai(x,z*)-fai(x,z^)
    U_x_zstar_zhat = fai_x_z_star{ind} - fai_x_z_hat;
    % sum( ��(x,z*,z^) )
    sum_U = U_x_zstar_zhat;
    % ����ÿһ����ÿ����������ʧ����
    sample_loss(t, ind) = delta_zstar_zhat;
    % =================================================================== %
    % sum( ��(z*,z^) ��
    sum_delta = sum(sample_loss(t,:));
    aver_loss(t) = sum_delta/1; 
    fprintf('      ��ǰ������ʧ������(z*,z^):\t%f\n', aver_loss(t));

    %% 3. ������Ų���������Wavg
    
    % ����ͷ��� lambda ��
    lambda = 1e-2;
	
    Ws = sum_U/(lambda*N);
    ls = sum_delta/N;
    
    % ����gap��gap��ֵ��������lambda����仯���޷�ȷ����������˻�����loss��gap�ȽϺ���
    gap_cur(t) = lambda*(Wi- Ws)'*W{t}- Li+ ls;
    % ���㲽��gamma
%     gamma(t) = gap_cur(t)/(lambda*norm(Wi- Ws)^2); % line search
    gamma(t) = 2*N/(2*N + t-1); % ��ͨ����
    
    % ���� wi��Li�������º��w������ W{t+1,ind}��
    Wi_old = Wi;
    Li_old = Li;
    Wi = (1- gamma(t))*Wi_old+ gamma(t)*Ws;
    Li = (1- gamma(t))*Li_old+ gamma(t)*ls;
    
    % ���� w��L�������º��w������ W{t+1,N+1}��
    W{t+1} = W{t} + Wi - Wi_old;
    L(t+1) = L(t) + Li - Li_old; % ����L�ƺ�ûʲô�ã�
    
    % ����Wavg
    Wavg{t+1} = (t-1)/(t+1)*Wavg{t} + 2/(t+1)*W{t+1};

    fprintf('      ��ż��϶gap:\t%f\n', gap_cur(t));
    % ==================================== %
    % ��¼ʱ��
    time(t) = toc;
%     fprintf('      ���Ƽ�϶ ��:\t%f\n', gap_cur);
    fprintf('      ʱ�仨��:\t%1.2f s\n', time(t)); 
    
    % �������ݣ���ֹ�ڴ治��ʱ���ݶ�ʧ
% 	save('E:\datasets\first editon\training datasets\N2DL-HeLa\01_2-16_track\�ṹ��ѧϰ\A_B_New.mat',...
% 	'N','s_frame','e_frame',... % 3��������ر���
% 	'iter','gap','gap_cur','A','B','W','J','t','time','sample_loss','aver_loss'); % 11��ȫ�ֱ���

end

% ѭ����ɣ���ӡ��Ϣ
if ls*N <= gap
    disp('  �ҵ��˵�ǰgap�µ����Ž⣬�㷨��ֹ');
    % �������ŷ��䷽����w
    t_best = t;
    w_best = Wavg{t_best}; % W{t}����ȡ����� gap ֵ���Ǹ�w
    gap_best = ls*N;      
else
    disp('  �ﵽ���ѭ���������㷨��ֹ');
    t_best = find(aver_loss==min(aver_loss(aver_loss~=0))); %  �ҵ���������ʧ��С���Ǹ�w��Ϊ w_best
    w_best = Wavg{t_best};
    gap_best = aver_loss(t_best);
end
% �������w�����ڲ�������֡����
% save([ trackpath, '\�ṹ��ѧϰ\SSVM_Best_W_New.mat'], 'w_best');

fprintf('\n\tt_best:\t%d\n', t_best);
fprintf('\tgap_best:\t%f\n', gap_best);
fprintf('\ttime consumption:\t%0.2f min\n', sum(time)/60);   
w_for_excel = w_best';

plot(aver_loss, '-*');
% �Եõ����������߽��б���
if 0
    name = 'loss_15_10_y';
    lossdir = [ trackpath, '\ѵ�������¼\BCFWavg_New\'];
    save([lossdir, name, '.mat'], 'aver_loss','sample_loss','w_best','Wavg');
    saveas(1, [lossdir, name, '.fig']);
end







