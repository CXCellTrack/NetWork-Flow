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
    w = rand(42,1);
%     w = zeros(42,1);
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

% ============================== ȫ�ֱ��� =============================== %
% ȫ���������ն��������������ֲ��������豣��

% ��ʼ���� ����A��B��ѭ���������� iter����϶��ֵ gap
iter = 50;
gap = 0.0010; % ���� O(1/gap) �������ٶȣ�Ӧ���ڰ�ѭ���������
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

tic;
% ���� fai(x,z)�͡�(z*,z)����������̱���
for ind=1:N
    disp('  ==========================');
    disp(['  Ԥ��������',num2str(ind),'��ѵ������...']);
    % ----------------------------------------- %
    % ������¼����̱�����Ԥ�ȼ���� fai(x,z)�� ��(z*,z)
    [ fij{ind} fit{ind} fid{ind} fiv{ind} fmj{ind} fsj{ind} phi_x_z{ind} sum_cost{ind} ] =...
        CXSL_Calculate_Fai_And_Loss( s_frame(ind), e_frame(ind) );
    % ����Լ������ F������ CXSL_Calculate_Constraint_New_Conflict �������
    % �� BundleMethod_Output_Test �е�ͬ������һ��
    % ----------------------------------------- %
    % 2015.7.6 ʹ�����µ�ì��Լ������22ì��Լ����
    [ F{ind} ] = CXSL_Calculate_Constraint_New_Conflict( 'training', true, s_frame(ind), e_frame(ind),...
        fij{ind}, fit{ind}, fid{ind}, fiv{ind}, fmj{ind}, fsj{ind} );
    % ----------------------------------------- %
	% �����׼���е�fai(x,z*)
	[ phi_x_z_star{ind} ] = CXSL_Calculate_fai_x_zstar_New( s_frame(ind), e_frame(ind), 'star'); 
    % ----------------------------------------- %
end
toc;

%% ����ǰѭ������tС�����ޣ���gap������Ҫ��ʱ������ѭ�����㣬�������󾫶Ȼ��������޸�gap��iter�������cell����

options = sdpsettings('verbose', 0, 'solver', 'cplex', 'saveduals', 0); % cplex���÷ŵ�ѭ����

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
    
    % ����ͷ��� lambda �ˣ����ԽС����wԽ�������ٶ��ƺ���ӿ죩
    lambda = 1e-6;
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
    gap_best = gap_cur(t);      
else
    disp('  �ﵽ���ѭ���������㷨��ֹ');
    t_best = find(aver_loss==min(aver_loss(aver_loss~=0))); %  �ҵ���������ʧ��С���Ǹ�w��Ϊ w_best
    w_best = W{t_best};
    gap_best = aver_loss(t_best);
end
% �������w�����ڲ�������֡����
save([ trackpath, '\�ṹ��ѧϰ\SSVM_Best_W_New.mat'], 'w_best');

% fprintf('\tw:\t%f\n', w_best);
fprintf('\tt_best:\t%d\n', t_best);
fprintf('\tgap_cur:\t%f\n', gap_best);
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








