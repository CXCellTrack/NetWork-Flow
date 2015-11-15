% ======================================================================= %
%
% 这个是 SSVM 训练的主函数 2015.6.17
% 按照 active structured learning 中的伪代码编写（Fig.4）
% 主要分为3个步骤：
%   1. 调用 CXSL_ILP 计算当前w下的最佳分配方案 z^
%
%   2. 计算梯度 U(x,z*,z^)=phi(x,z*)-phi(x,z^)，并得到a与b的值（公式12、13）
%
%   3. 通过a、b解方程14求出更新后的 w和 kexi，计算gap的大小
%
%   运行中发现有内存不足的问题，主要是没有连续的内存，因此通过预分配变量空间、
% 每运行50次 clear 一下变量并重新载入来解决
% 
%
% ======================================================================= %

% 这个脚本是完全按照 BCFW 论文中的算法4编写
clear;close all;

[ ~, trackpath ] = getpath( 'training' );
% 载入 CXSL_Test_Linear_all 中计算好的 w 作为初始值（实际发现效果并不好）
load([ trackpath, '\结构化学习\initial_w_New.mat']);
% 注意 w 的顺序不能乱
w = [ wij, wit, wid, wiv, wmj, wsj ]';
clear wij wit wid wiv wmj wsj
% 也可以选择随机的w或全0的w
if 0
    w = zeros(numel(w), 1);
end

% 定义样本个数 N 和 单个样本中的帧数 frame
N = 8;
frame = 10;
s_frame = zeros(N,1);
e_frame = zeros(N,1);
% 目前有gt的帧数，对随机取样有影响
gt_frame = 80;

% 选择取样方式
% 1为接龙取样，2为滑窗取样，3为随机取样
sample_method = 1;
switch sample_method
    case 1
        % 取样方法1：接龙取样
        % 按 1-5, 6-10, 11-15, 16-20 这样的方法取样本
        for ind=1:N
            s_frame(ind) = (ind - 1)*frame + 1;
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
iter = 200;
gap = 0.0010; % 按照 O(1/gap) 的收敛速度，应该在百循环左右完成

%% ------------------ 线性w所用到的变量 ---------------------- %
W = cell(iter,1); % W 存放综合权值w
Wi = cell(iter,N); % Wi存放样本权值w
W{1} = w; % 全体样本的W需要设定初值w
for i=1:N
    Wi{1,i} = w; % Wi 存放每次循环中特定样本更新后的Wi
end

L = zeros(iter,1); % L 存放综合L
Li = zeros(iter,N);
for i=1:N
    Li(1,i) = 0;
end
ls = 1; % 样本平均损失函数

%% ------------------ 循环中用到的变量 -------------------- %
gap_cur = zeros(iter,1); % 记录每次得到的gap
gamma = zeros(iter,1); % 步长gamma
t = 0;
time = zeros(iter,1); % 记录每次循环所用的时间

sample_loss = zeros(iter,N); % 记录每一轮中每个样本的损失函数
aver_loss = zeros(iter,1); % 记录每一轮中样本损失函数均值
% ======================================================================= %
% BCFW模式下每次只训练一个样本，因此这些原本需要累加的量不需要了
% phi_x_z_hat = cell(N,1);
% delta_zstar_zhat = zeros(N,1);
% U_x_zstar_zhat = cell(N,1);
% ======================================================================== %
% 循环求解部分参数设置
options = sdpsettings('verbose', 0, 'solver', 'cplex'); % cplex设置放到循环外
rng(0); % 含有随机选择部分，需要设定种子
random = 0; % 变量random作为一个flag，为1时是随机抽样，为0时是按顺序抽样
ind = 0;

disp('  预计算目标函数和约束条件...');

%% 提前计算好 phi(x,z^) phi(x,z*)和△(z*,z^)，还有约束条件，循环中组装目标函数，再求解
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
% 计算 phi(x,z)和△(z*,z)，分配好流程变量
for ind=1:N
    disp('  ==========================');
    disp(['  预计算样本',num2str(ind),'的训练数据...']);
    % ----------------------------------------- %
    % 分配各事件流程变量，预先计算好 phi(x,z)和 △(z*,z)
    [ fij{ind} fit{ind} fid{ind} fiv{ind} fmj{ind} fsj{ind} phi_x_z{ind} sum_cost{ind} ] =...
        CXSL_Calculate_phi_And_Loss( s_frame(ind), e_frame(ind) );
    % 计算约束条件 F，调用 CXSL_Calculate_Constraint_New_Conflict 这个函数
    % 与 BundleMethod_Output_Test 中的同名函数一样
    % ----------------------------------------- %
    % 2015.7.6 使用了新的矛盾约束规则（22矛盾约束）
    [ F{ind} ] = CXSL_Calculate_Constraint_New_Conflict( 'training', true, s_frame(ind), e_frame(ind),...
        fij{ind}, fit{ind}, fid{ind}, fiv{ind}, fmj{ind}, fsj{ind} );
    % ----------------------------------------- %
	% 计算标准答案中的phi(x,z*)
	[ phi_x_z_star{ind} ] = CXSL_Calculate_phi_x_zstar_New( s_frame(ind), e_frame(ind), 'star'); 
    % ----------------------------------------- %
end
toc;

%% 当当前循环次数t小于上限，且gap不符合要求时，进行循环计算，若想增大精度或轮数，修改gap和iter再运行此cell即可

while t < iter && ls*N >= gap
    t = t + 1;
    
    % 记录下每次循环所用的时间
    t_start = clock;
    disp('  ==========================');
    disp(['  开始第 ', num2str(t), ' 轮循环...']);
    
    % 可以选择随机抽一个样本或是按顺序来
    if random
        ind = randi(N); % 选中第ind个样本作为训练样本
    else
        ind = ind + 1;
        if ind==N + 1
            ind = 1;
        end
    end
    
    %% 1. 计算给定 w 下，当前选定样本的最佳分配（方程（10））
    
    disp(['      计算样本 ', num2str(ind), '...']); 
    % 临时组建目标函数并求解
    object_function = dot(W{t}, phi_x_z{ind}) + sum_cost{ind};
    sol = solvesdp( F{ind}, -object_function, options );

    % 输出得到的各个变量的值
    if sol.problem == 0      
        phi_x_z_hat = value(phi_x_z{ind});
        delta_zstar_zhat = value(sum_cost{ind});
    else
        sol.info
        yalmiperror(sol.problem)
    end
   
    disp('      更新权向量 w...');
    
    %% 2. 计算 ψ(x,z*,z^)=phi(x,z*)-phi(x,z^) 梯度
    % 用 U 来代表这个符号 ψ
    
%     for ind=1:N
        % 梯度 论文中phi(x,z*)-phi(x,z^)
        U_x_zstar_zhat = phi_x_z_star{ind} - phi_x_z_hat;
        % sum( ψ(x,z*,z^) )
        sum_U = U_x_zstar_zhat;
        % 保存每一轮中每个样本的损失函数
        sample_loss(t, ind) = delta_zstar_zhat;
%     end
    % =================================================================== %
    % sum( △(z*,z^) ）
    sum_delta = sum(sample_loss(t,:));
    aver_loss(t) = sum_delta; 
    fprintf('      当前样本损失函数△(z*,z^):\t%f\n', aver_loss(t));

    %% 3. 求解最优步长来更新w   
    % 定义惩罚项 lambda λ
    lambda = 1e-2;
    Ws = sum_U/(lambda*N);
    ls = sum_delta/N;  
    % 计算gap，gap的值随样本、lambda都会变化，无法确定下来，因此还是用loss做gap比较合适
    % 计算步长gamma
    linesearch = 1;
    if linesearch
        gap_cur(t) = lambda*(Wi{t,ind}- Ws)'*W{t}- Li(t,ind)+ ls;
        gamma(t) = gap_cur(t)/(lambda*norm(Wi{t,ind}- Ws)^2);
        gamma(t) = max([0, min([gamma(t),1])]);
    else
        gamma(t) = 2*N/(2*N + t-1);
    end
    % 更新 wi和Li，将更新后的w保存在 W{t+1,ind}中
    Wi{t+1,ind} = (1- gamma(t))*Wi{t,ind} + gamma(t)*Ws;
    Li(t+1,ind) = (1- gamma(t))*Li(t,ind) + gamma(t)*ls;
    % 对于此轮没轮到的样本，将其Wi和Li带入下一轮中
    sample_not_used = mysetdiff(1:N,ind);
    for jj=sample_not_used
        % jj为没被用到的样本编号
        Wi{t+1,jj} = Wi{t,jj};
        Li(t+1,jj) = Li(t,jj); % 直接将其Wi带入下一轮
    end
    
    % 更新 w和L，将更新后的w保存在 W{t+1,N+1}中
    W{t+1} = W{t} + Wi{t+1,ind} - Wi{t,ind};
    L(t+1) = L(t) + Li(t+1,ind) - Li(t,ind); % 计算L似乎没什么用？

    fprintf('      对偶间隙gap:\t%f\n', gap_cur(t));
    % ==================================== %
    % 记录时间
    time(t) = etime(clock, t_start);
%     fprintf('      近似间隙 ε:\t%f\n', gap_cur);
    fprintf('      时间花费:\t%1.2f s\n', time(t)); 

end

% 循环完成，打印信息
if ls*N <= gap
    disp('  找到了当前gap下的最优解，算法终止');
    
    % 保存最优分配方案和w
    t_best = t;
    w_best = W{t_best}; % W{t-1}才是取到最佳 gap 值的那个w，随后更新得到的 W{t} gap可能增大了
    gap_best = ls*N;      
else
    disp('  达到最大循环次数，算法终止');
    t_best = find(aver_loss==min(aver_loss(aver_loss~=0))); %  找到过程中损失最小的那个w作为 w_best
    w_best = W{t_best};
    gap_best = aver_loss(t_best);
end
% 保存最佳w，用于测试其他帧精度
save([ trackpath, '\结构化学习\SSVM_Best_W_New.mat'], 'w_best');

fprintf('\n\tt_best:\t%d\n', t_best);
fprintf('\tgap_best:\t%f\n', gap_best);
fprintf('\ttime consumption:\t%0.2f min\n', sum(time)/60);   
w_for_excel = w_best';

plot(aver_loss, '-*');
% 对得到的收敛曲线进行保存
if 0
    name = 'loss_8_10_cons1235';
    lossdir = [ trackpath, '\训练结果记录\BCFW_paper\'];
    mkdir(lossdir);
    save([lossdir, name, '.mat'], 'w','linesearch','sample_loss','lambda','w_best','W');
    saveas(1, [lossdir, name, '.fig']);
end




