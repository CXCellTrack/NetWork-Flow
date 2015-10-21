% ======================================================================= %
%
% 这个是 SSVM 训练的主函数 2015.6.17
% 按照 active structured learning 中的伪代码编写（Fig.4）
% 主要分为3个步骤：
%   1. 调用 CXSL_ILP 计算当前w下的最佳分配方案 z^
%
%   2. 计算梯度 U(x,z*,z^)=fai(x,z*)-fai(x,z^)，并得到a与b的值（公式12、13）
%
%   3. 通过a、b解方程14求出更新后的 w和 kexi，计算gap的大小
%
%   运行中发现有内存不足的问题，主要是没有连续的内存，因此通过预分配变量空间、
% 每运行50次 clear 一下变量并重新载入来解决
% 
%
% ======================================================================= %
clear;close all;
% 载入 CXSL_Test_Linear_all 中计算好的 w 作为初始值（实际发现效果并不好）
if 1
    load('E:\datasets\first editon\training datasets\N2DL-HeLa\01_2-16_track\结构化学习\initial_w_New.mat');
    % 注意 w 的顺序不能乱
    w = [ wij, wit, wid, wiv, wmj, wsj ]';
    clear wij wit wid wiv wmj wsj
else
    % 随机选取w，种子控制 w 可复现
    rng(0); 
    w = rand(42,1);
%     w = zeros(42,1);
end
% 定义样本个数 N 和 单个样本中的帧数 frame
N = 1;
frame = 3;
s_frame = zeros(N,1);
e_frame = zeros(N,1);
% 目前有gt的帧数，对随机取样有影响
gt_frame = 70;

% 选择取样方式
% 1为接龙取样，2为滑窗取样，3为随机取样
sample_method = 2;
switch sample_method
    case 1
        % 取样方法1：接龙取样
        % 按 1-5, 5-9, 9-13, 13-17 这样的方法取样本
        for ind=1:N
            s_frame(ind) = (ind - 1)*(frame - 1) + 1;
            e_frame(ind) = s_frame(ind) + frame - 1;
        end

    case 2
        % 取样方法2：滑窗取样
        % 按 1-5，2-6，3-7 这样的方法取样本
        for ind=1:N
            s_frame(ind) = ind;
            e_frame(ind) = s_frame(ind) + frame - 1;
        end
        
    case 3
        % 取样方法3：随机取样
        rng(0);
        s_frame = randi([1 gt_frame-frame+1], [N 1]);
        e_frame = s_frame + frame - 1;        
end

% 打印样本信息
disp(['  训练共选取', num2str(N), '个样本：']);
for ind=1:N
    disp(['  ', num2str(s_frame(ind)), '——', num2str(e_frame(ind)), '帧...']);
end

% ============================== 全局变量 =============================== %
% 全部变量最终都被保存下来，局部变量无需保存

% 初始化： 定义A、B、循环次数上限 iter、间隙阈值 gap
iter = 300;
gap = 0.05; % 按照 O(1/gap) 的收敛速度，应该在百循环左右完成
gap_cur = zeros(iter,1); % 记录每次得到的gap
gap_cur(1) = gap* 2;

time = zeros(iter,1); % 记录每次循环所用的时间
sample_loss = zeros(iter,N); % 记录每一轮中每个样本的损失函数
now_loss = zeros(iter,1); % 记录每一轮中样本损失函数均值
gamma = zeros(iter,1); % 步长gamma

W = cell(iter,1); % W 存放综合权值w
Wi = cell(iter,N); % Wi 存放每次循环中特定样本更新后的Wi
Wi_old = cell(1,N); % Wi_old 存放每次循环中特定样本更新前的Wi
Ws = cell(iter,N);
L = zeros(iter,1); % L 存放综合L
Li = zeros(iter,N); % Li 存放每次选中的样本更新后的L
Li_old = zeros(1,N); % Li_old 存放每次选中的样本更新前的L
ls = 1; % 样本平均损失函数

W{1} = w; % 全体样本的W需要设定初值w
for i=1:N
    Wi{1,i} = w; % 单个样本的Wi也设定初值为w，否则会导致单样本BCFW与FW之间的结果不一致
    Wi_old{i} = w;
end

% ------------------- 以下为MP-BCFW增加的变量 ----------------- %
N_set = 1000; % 工作集大小（设为一个大数）
M_pass = 1000; % 运行近似oracle的最多次数，实际为动态决定，当达到最大后则强行停止
T_last = 10; % 最后T次循环，用于移除Wi
Fobj = zeros(iter,N); % 存放每个样本的目标函数

t = 1;
% ======================================================================= %

fai_x_z_hat = cell(N,1);
delta_zstar_zhat = zeros(N,1);

% ============ 步骤2的变量 ========== %
U_x_zstar_zhat = cell(N,1);
% ======================================================================== %

disp('  预计算目标函数和约束条件...');
%% 提前计算好 fai(x,z^) fai(x,z*)和△(z*,z^)，还有约束条件，循环中组装目标函数，再求解
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
% 计算 fai(x,z)和△(z*,z)，分配好流程变量
for ind=1:N
    % ----------------------------------------- %
    % 分配各事件流程变量，预先计算好 fai(x,z)和 △(z*,z)
    [ fij{ind} fit{ind} fid{ind} fiv{ind} fmj{ind} fsj{ind} fai_x_z{ind} sum_cost{ind} ] =...
        CXSL_Calculate_Fai_With_Loss_New( s_frame(ind), e_frame(ind) );
    % 计算约束条件 F，调用 CXSL_Calculate_Constraint_New_Conflict 这个函数
    % 与 BundleMethod_Output_Test 中的同名函数一样
    % ----------------------------------------- %
    % 2015.7.6 使用了新的矛盾约束规则（22矛盾约束）
    [ F{ind} ] = CXSL_Calculate_Constraint_New_Conflict( s_frame(ind), e_frame(ind),...
        fij{ind}, fit{ind}, fid{ind}, fiv{ind}, fmj{ind}, fsj{ind} );
    % ----------------------------------------- %
	% 计算标准答案中的fai(x,z*)
	[ fai_x_z_star{ind} ] = CXSL_Calculate_fai_x_zstar_New( s_frame(ind), e_frame(ind), 'star'); 
    % ----------------------------------------- %
end
toc;

%% 当当前循环次数t小于上限，且gap不符合要求时，进行循环计算，若想增大精度或轮数，修改gap和iter再运行此cell即可

options = sdpsettings('verbose', 0, 'solver', 'cplex', 'saveduals', 0); % cplex设置放到循环外
rng(0); % 含有随机选择部分，需要设定种子

while t <= iter && ls*N >= gap

    % 记录下每次循环所用的时间
    t_start = clock;
    disp('  ==========================');
    disp(['  开始第 ', num2str(t), ' 轮循环...']);
    
	%% step1. 更新一次精确oracle
    tic;
    
    	%% 1. 计算给定w下的精确oracle
    
    ind = randi(N);
    disp(['      计算样本 ', num2str(ind), '...']); 
    % 临时组建目标函数并求解
    object_function = dot(W{t}, fai_x_z{ind}) + sum_cost{ind};
    sol = solvesdp( F{ind}, -object_function, options );

    % 输出得到的各个变量的值
    if sol.problem == 0      
        fai_x_z_hat{ind} = value(fai_x_z{ind});
        delta_zstar_zhat(ind) = value(sum_cost{ind});
        Fobj(t,ind) = value(object_function); % 这一次计算oracle得到的目标函数值
    else
        sol.info
        yalmiperror(sol.problem)
    end
   
    disp('      更新权向量 w...');
    
        %% 2. 计算 ψ(x,z*,z^)=fai(x,z*)-fai(x,z^) 梯度
    % 用 U 来代表这个符号 ψ
    
%     for ind=1:N
        % 梯度 论文中fai(x,z*)-fai(x,z^)
        U_x_zstar_zhat{ind} = fai_x_z_star{ind} - fai_x_z_hat{ind};
        % sum( ψ(x,z*,z^) )
        sum_U = U_x_zstar_zhat{ind};
        % 保存每一轮中每个样本的损失函数
        sample_loss(t, ind) = delta_zstar_zhat(ind);
%     end
    % =================================================================== %
    % sum( △(z*,z^) ）
    sum_delta = sum(sample_loss(t,:));
    now_loss(t) = sum_delta; 
    fprintf('      当前样本损失函数△(z*,z^):\t%f\n', sum_delta);

        %% 3. 求解最优步长来更新w
    
    % 定义惩罚项 lambda λ
    lambda = 1e-2;
	
    Ws{t,ind} = sum_U/(lambda*N);
    ls = sum_delta/N;
    
    % 计算gap，gap的值随样本、lambda都会变化，无法确定下来，因此还是用loss做gap比较合适
    gap_cur(t+1) = lambda*(Wi_old{ind}- Ws{t,ind})'*W{t}- Li_old(ind)+ ls;
    % 计算步长gamma
    gamma(t) = gap_cur(t+1)/(lambda*norm(Wi_old{ind} - Ws{t,ind})^2);
    
    % 更新 wi和Li，将更新后的w保存在 W{t+1,ind}中
    Wi{t+1,ind} = (1- gamma(t))*Wi_old{ind}+ gamma(t)*Ws{t,ind};
    Li(t+1,ind) = (1- gamma(t))*Li_old(ind)+ gamma(t)*ls;

    % 更新 w和L，将更新后的w保存在 W{t+1,N+1}中
    W{t+1} = W{t} + Wi{t+1,ind} - Wi_old{ind};
    L(t+1) = L(t) + Li(ind) - Li_old(ind); % 计算L似乎没什么用？
    Wi_old{ind} = Wi{t+1,ind};
    Li_old(ind) = Li(t+1,ind);
    
    fprintf('      对偶间隙gap:\t%f\n', gap_cur(t+1));

        %% 4. 判断工作集长度是否超过设定值Nset
        N_work = 0;
        for uu=1:t
            if ~isempty(Wi{uu,ind})
                N_work = N_work + 1;
            end
        end
        % 超出工作集长度时需要删去最不活跃的一个Wi
        if N_work>N_set
            %
            % 未完待续
            %
        end
      
    %% step2. 求解若干次近似oracle
    t_exact = toc;
    
        %% 1. 在工作集中找到能使目标函数值最大的那个W
        tbest = find(Fobj==max(Fobj));
        % 更新w
        

        
        
        
    % 进入下一轮循环
    t = t + 1;  
    % ==================================== %
    % 记录时间
    t_end  = clock;
    time(t) = etime(t_end, t_start);
    fprintf('      时间花费:\t%1.2f s\n', time(t)); 
end

% 循环完成，打印信息
if ls*N <= gap
    disp('  找到了当前gap下的最优解，算法终止');
    
    % 保存最优分配方案和w
    t_best = t - 1;
    w_best = W{t-1}; % W{t-1}才是取到最佳 gap 值的那个w，随后更新得到的 W{t} gap可能增大了
    gap_best = ls*N;      
else
    disp('  达到最大循环次数，算法终止');
    t = t - 1;
    t_best = find(now_loss==min(now_loss)); %  找到过程中损失最小的那个w作为 w_best
    w_best = W{t_best};
    gap_best = now_loss(t_best);
end
% 保存最佳w，用于测试其他帧精度
save('E:\datasets\first editon\training datasets\N2DL-HeLa\01_2-16_track\结构化学习\SSVM_Best_W_New.mat', 'w_best');

figure;plot(now_loss);
% subplot(211);  %aver_loss = aver_loss(1:t); 
% subplot(212); plot(gap_cur); %gap_cur = gap_cur(1:t); 
% fprintf('\tw:\t%f\n', w_best);
fprintf('\n\tt_best:\t%d\n', t_best);
fprintf('\tgap_best:\t%f\n', gap_best);
fprintf('\ttime consumption:\t%0.2f min\n', sum(time)/60);   
% w_for_excel = w_best';









