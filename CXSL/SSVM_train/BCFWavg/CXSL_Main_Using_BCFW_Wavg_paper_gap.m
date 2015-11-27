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

% 此方法每个样本都保存一个Wi，使得样本更新的时候总是从上一次的Wi=0开始
clear;close all;

[ ~, trackpath ] = getpath( 'training' );
% 载入 CXSL_Test_Linear_all 中计算好的 w 作为初始值（实际发现效果并不好）
load([ trackpath, '\结构化学习\initial_w_New.mat']);
% 注意 w 的顺序不能乱
w = [ wij,bij, wit,bit, wid,bid, wiv,biv, wmj,bmj, wsj,bsj ]'; % 原版b
if 0
    % 第一个数据集
%     w = [-24.05027785	0.965959752	-0.209700235	0.023655542	-0.901444678	0.915527485	-0.723055368	0.78127324	-22.34659216	-3.45491283	-1.682414322	-5.355960441	-2.391659001	2.862181421	-7.382944338	8.382838223	1.94377663	-0.451290137	-1.07738777	-4.844423375	-1.122913059	-0.801496889	3.907101647	-11.61160994	3.710115534	0.998335816	4.252699702	0.790594494	1.207125853	3.799458373	1.390618031	5.18991389	1.129864864	0.673380786	-2.076937813	-1.97433464	-1.980221778	-0.051210814	0.597328997	-3.897482158]';
    % 第二个数据集
%     w = [-18.21315239	2.048551055	0.090611096	-0.07830978	-0.681768441	0.091705287	-0.284558766	0.113666465	-16.77170209	-2.820207584	-1.606735489	-2.170929556	-2.000511632	2.42450433	-4.406444861	10.34098417	1.758814312	-0.819906672	-2.159095585	-4.969572233	1.29607646	0.202318113	3.379177651	-14.25716614	7.109097539	3.366559674	2.10659084	-2.499899814	-3.9849466	2.501397849	1.351917771	4.699879458	-0.525793863	-0.261628668	-4.945586679	-1.846739083	-4.998199696	-0.003648138	1.737582541	-8.35761154]';
    % 第三个数据集
    w = [-17.69956928	1.61055833	-0.001787913	0.199592206	-0.287005286	-1.058810184	-0.194715268	0.170821175	-15.86943474	-2.459629861	-1.102503407	-1.233255838	-1.209682446	-2.849135318	-5.932152347	3.635597276	3.112043952	0.905772748	-0.112320949	-0.658275898	4.373948249	0.504445114	-0.211398811	-3.571706443	0	0	0	0	0	0	0	0	0	0	-2.873279295	-0.074622812	-1.438747225	-0.381464188	-2.838608426	-5.896985464]';
else
    % initial_w是增广的
    w = zeros(numel(w), 1);
end

% 定义样本个数 N 和 单个样本中的帧数 frame
N = 5;
frame = 13;
s_frame = zeros(N,1);
e_frame = zeros(N,1);
% 目前有gt的帧数，对随机取样有影响
gt_frame = 65;

% 选择取样方式
sample_method = 2;
switch sample_method
    case 1 % 取样方法1：接龙取样（原方法）
        % 按 1-5, 6-10, 12-15, 16-20 这样的方法取样本
        for ii=1:N
            s_frame(ii) = (ii - 1)*frame + 1;
            e_frame(ii) = s_frame(ii) + frame - 1;
        end

    case 2 % 取样方法2：重叠接龙取样（new sample method）
        % 按 1-5, 5-10, 10-15, 15-20 这样的方法取样本
        for ii=1:N
            s_frame(ii) = (ii - 1)*frame;
            e_frame(ii) = s_frame(ii) + frame;
        end
        s_frame(1) = 1;
        
    case 3 % 取样方法3：滑窗取样
        % 按 1-5，2-6，3-7 这样的方法取样本
        for ind=1:N
            s_frame(ind) = ind;
            e_frame(ind) = s_frame(ind) + frame - 1;
        end
       
    case 4 % 取样方法4：随机取样
        rng(0);
        s_frame = randi([1 gt_frame-frame+1], [N 1]);
        e_frame = s_frame + frame - 1;        
end

% 打印样本信息
disp(['  训练共选取', num2str(N), '个样本：']);
for ind=1:N
    disp(['  ', num2str(s_frame(ind)), '——', num2str(e_frame(ind)), '帧...']);
end

% 全部变量最终都被保存下来，局部变量无需保存

% 初始化： 定义A、B、循环次数上限 iter、间隙阈值 gap
iter = 1000;
gap = 0.0010; % 按照 O(1/gap) 的收敛速度，应该在百循环左右完成

%% ============================== 全局变量 =============================== %
gap_cur = zeros(iter,1); % 记录每次得到的gap
gamma = zeros(iter,1); % 步长gamma

W = cell(iter,1); % W 存放综合权值w
Wavg = cell(iter,1); % W 存放综合权值w
Wi = cell(iter,N); % Wi存放样本权值w
Wavg{1} = w;
W{1} = w; % 全体样本的W需要设定初值w
for i=1:N
    Wi{1,i} = w; % Wi 存放每次循环中特定样本更新后的Wi
end

L = zeros(iter,1); % L 存放综合L
Li = zeros(iter,N);
for i=1:N
    Li(1,i) = 0;
end

t = 0;
time = zeros(iter,1); % 记录每次循环所用的时间
sample_loss = zeros(iter,N); % 记录每一轮中每个样本的损失函数
aver_loss = zeros(iter,1); % 记录每一轮中样本损失函数均值
% ======================================================================= %
phi_x_z_hat = cell(N,1);
delta_zstar_zhat = zeros(N,1);
% U_x_zstar_zhat = cell(N,1);
% ======================================================================== %
% 循环求解部分参数设置
options = sdpsettings('verbose', 0, 'solver', 'gurobi','cachesolvers',1); % cplex设置放到循环外
rng(0); % 含有随机选择部分，需要设定种子
random = 1; % 变量random作为一个flag，为1时是随机抽样，为0时是按顺序抽样
ind = 0;

disp('  预计算目标函数和约束条件...');

%% 提前计算好 phi(x,z^) phi(x,z*)和△(z*,z^)，还有约束条件，循环中组装目标函数，再求解
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
% 计算 phi(x,z)和△(z*,z)，分配好流程变量
for ii=1:N
    disp('  ==========================');
    disp(['  预计算样本',num2str(ii),'的训练数据...']);
    % ----------------------------------------- %
    % 分配各事件流程变量，预先计算好 phi(x,z)和 △(z*,z)
    [ fij{ii} fit{ii} fid{ii} fiv{ii} fmj{ii} fsj{ii} phi_x_z{ii} sum_cost{ii} sum_cost_all{ii} ] =...
        CXSL_Calculate_phi_And_Loss( w, s_frame(ii), e_frame(ii) );
    % 计算约束条件 F，调用 CXSL_Calculate_Constraint_New_Conflict 这个函数
    % 与 BundleMethod_Output_Test 中的同名函数一样
    % ----------------------------------------- %
    % 2015.7.6 使用了新的矛盾约束规则（22矛盾约束）
    [ F{ii} ] = CXSL_Calculate_Constraint_New_Conflict( 'training', use_op_cons, s_frame(ii), e_frame(ii),...
        fij{ii}, fit{ii}, fid{ii}, fiv{ii}, fmj{ii}, fsj{ii} );
    % ----------------------------------------- %
	% 计算标准答案中的phi(x,z*)
	[ phi_x_z_star{ii} ] = CXSL_Calculate_phi_x_zstar_New( w, s_frame(ii), e_frame(ii), 'star'); 
    % ----------------------------------------- %
end
toc;

%% 当当前循环次数t小于上限，且gap不符合要求时，进行循环计算，若想增大精度或轮数，修改gap和iter再运行此cell即可
% 定义惩罚项 lambda λ
n4gap = 2; % 每N*n4gap轮循环后计算下gap（n4gap>1）
lambda = 1e-2;
linesearch = 0;
usecostall = 0;
if usecostall
    disp('当前选择的损失中包含了虚景！');
    sum_cost = sum_cost_all;
end

while t < iter % && ls*N >= gap
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
    if mod(t,n4gap*N)~=0 && t~=1
        object_function = dot(Wavg{t}, phi_x_z{ind}) + sum_cost{ind};
        sol = solvesdp( F{ind}, -object_function, options );

        % 输出得到的各个变量的值
        if sol.problem == 0      
            phi_x_z_hat{ind} = value(phi_x_z{ind});
            delta_zstar_zhat(ind) = value(sum_cost{ind});
        else
            sol.info
            yalmiperror(sol.problem)
        end
    else
        disp('      计算一次全体样本的损失...');
        for ii=1:N
            disp(['        样本',num2str(ii),'...']);
            object_function = dot(Wavg{t}, phi_x_z{ii}) + sum_cost{ii};
            sol = solvesdp( F{ii}, -object_function, options );
            % 输出得到的各个变量的值
            if sol.problem == 0      
                phi_x_z_hat{ii} = value(phi_x_z{ii});
                delta_zstar_zhat(ii) = value(sum_cost{ii});
            else
                sol.info
                yalmiperror(sol.problem)
            end
        end
    end
        
    disp('      更新权向量 w...');
    
    %% 2. 计算 ψ(x,z*,z^)=phi(x,z*)-phi(x,z^) 梯度
    % 用 U 来代表这个符号 ψ
    sum_U = 0;
    if mod(t,n4gap*N)~=0 && t~=1
        U_x_zstar_zhat = phi_x_z_star{ind} - phi_x_z_hat{ind};
        sum_U = U_x_zstar_zhat;
        % 保存每一轮中每个样本的损失函数
        sample_loss(t, ind) = delta_zstar_zhat(ind);
    else
        % 计算所有样本的平均损失
        for ii=1:N
            % 梯度 论文中phi(x,z*)-phi(x,z^)
            U_x_zstar_zhat = phi_x_z_star{ii} - phi_x_z_hat{ii};
            % sum( ψ(x,z*,z^) )
            sum_U = sum_U + U_x_zstar_zhat;
            % 保存每一轮中每个样本的损失函数
            sample_loss(t, ii) = delta_zstar_zhat(ii);
        end
        
    end
    % =================================================================== %
    % sum( △(z*,z^) ）
    sum_delta = sum(sample_loss(t,:));
    aver_loss(t) = sum_delta;
    fprintf('      当前样本损失函数△(z*,z^):\t%f\n', sample_loss(t,ind));
    if mod(t,n4gap*N)==0 || t==1
        aver_loss(t) = sum_delta/N;
        fprintf('      当前pass的平均损失函数△(z*,z^):\t%f\n', aver_loss(t));
    end

    %% 3. 求解最优步长来更新Wavg
    Ws = sum_U/(lambda*N);
    ls = sum_delta/N;  
    % 计算gap，gap的值随样本、lambda都会变化，无法确定下来，因此还是用loss做gap比较合适
    if mod(t,n4gap*N)==0 || t==1
        gap_cur(t) = lambda*(W{t}- Ws)'*W{t}- L(t)+ ls;
    end
    % 计算步长gamma
    if linesearch
        tmp = lambda*(Wi{t,ind}- Ws)'*W{t}- Li(t,ind)+ ls;
        gamma(t) = tmp/(lambda*norm(Wi{t,ind}- Ws)^2);
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
    L(t+1) = L(t) + Li(t+1,ind) - Li(t,ind); % 计算L似乎没什么用？（有用！2015.11.8）
    
    % 更新Wavg
    Wavg{t+1} = (t-1)/(t+1)*Wavg{t} + 2/(t+1)*W{t+1};

    fprintf('      对偶间隙gap:\t%f\n', gap_cur(t));
    % ==================================== %
    % 记录时间
    time(t) = etime(clock, t_start);
    fprintf('      时间花费:\t%1.2f s\n', time(t)); 
    
end

% 循环完成，打印信息
if ls*N <= gap
    disp('  找到了当前gap下的最优解，算法终止');
    % 保存最优分配方案和w
    t_best = t;
    w_best = Wavg{t}; % W{t}才是取到最佳 gap 值的那个w
    loss_best = ls*N;      
else
    disp('  达到最大循环次数，算法终止');
    
    t_best = find(aver_loss==min(aver_loss(N+1:end))) %  找到过程中损失最小的那个w作为 w_best
    t_best = t_best(t_best>N);
    t_best = t_best(end);
    % gap小的学习损失不一定小，test误差也不一定小（还是依靠学习损失比较准）
%     t_best = find(gap_cur==min(gap_cur(gap_cur~=0)));
    w_best = Wavg{t_best};
    loss_best = aver_loss(t_best);
end
% 保存最佳w，用于测试其他帧精度
save([ trackpath, '\结构化学习\SSVM_Best_W_New.mat'], 'w_best');

fprintf('\n\tt_best:\t%d\n', t_best);
fprintf('\tloss_best:\t%f\n', loss_best);
fprintf('\ttime consumption:\t%0.2f min\n', sum(time)/60);   
w_for_excel = w_best';

% 绘制loss和gap
subplot(211)
    plot(aver_loss, '-*');hold on;
    gapt = 0:n4gap*N:iter;gapt(1) = 1;
    gaploss = aver_loss(gapt);
    plot(gapt,gaploss, 'r-o','linewidth',2);hold off;
subplot(212)
    gapc = gap_cur(gapt);
    plot(gapt,gapc, 'r-o','linewidth',2);
% 对得到的收敛曲线进行保存
if 0
    name = 'loss_5_13_cons35_cost1_init0p_noline_rng';
    lossdir = [ trackpath, '\训练结果记录\BCFWavg_paper\new_sample_method\'];
    mkdir(lossdir);
    save([lossdir, name, '.mat'], 'time','w','linesearch','use_op_cons','sample_loss','lambda','w_best','Wavg');
    saveas(1, [lossdir, name, '.fig']);
end







