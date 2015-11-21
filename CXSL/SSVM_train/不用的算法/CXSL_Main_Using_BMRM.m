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

[ ~, trackpath ] = getpath( 'training' );
load([ trackpath, '\结构化学习\initial_w_New.mat']);
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
sample_method = 1;
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

% 初始化： 定义A、B、循环次数上限 iter、间隙阈值 gap
iter = 50;
gap = 0.0010; % 按照 O(1/gap) 的收敛速度，应该在百循环左右完成

%% ============================== 全局变量 =============================== %
% 全部变量最终都被保存下来，局部变量无需保存
gap_cur = zeros(iter,1); % 记录每次得到的gap

A = cell(iter,1);
B = cell(iter,1);
W = cell(iter,1); % W 存放权值w
W{1} = w;
ls = 1;
t = 0;
time = zeros(iter,1); % 记录每次循环所用的时间

sample_loss = zeros(iter,N); % 记录每一轮中每个样本的损失函数
aver_loss = zeros(iter,1); % 记录每一轮中样本损失函数均值
% ======================================================================= %
% 增量学习：中间断开的话可以从数据中载入 A、B、w 等继续循环
% ============================== 局部变量 ================================ %
% 即每一次循环中这些变量的值都会被替换掉，预分配空间可以加快速度
% ============ 步骤1的变量 ========== %
phi_x_z_hat = cell(N,1);
delta_zstar_zhat = zeros(N,1);
psi_zstar_zhat = cell(N,1);
% ======================================================================= %
disp('  预计算目标函数和约束条件...');

%% 提前计算好 fai(x,z^) fai(x,z*)和△(z*,z^)，还有约束条件，循环中组装目标函数，再求解
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

%% 当当前循环次数t小于上限，且gap不符合要求时，进行循环计算，若想增大精度或轮数，修改gap和iter再允许此cell即可
options = sdpsettings('verbose', 0, 'solver', 'gurobi','cachesolvers',1); % cplex设置放到循环外
% 定义惩罚项 lambda λ，这个越小，则w越大（收敛速度似乎会加快）
lambda = 1e-8;
usecostall = 0;
if usecostall
    sum_cost = sum_cost_all;
end

while t < iter && ls >= gap

    t = t + 1;
    % 记录下每次循环所用的时间
    tic;
    disp('  ==========================');
    disp(['  开始第 ', num2str(t), ' 轮循环...']);
    
    %% 1. 计算给定 w 下，每个样本的最佳分配（方程（10））
    
    for ind=1:N
        disp(['      计算样本 ', num2str(ind), '...']); 
        % 临时组建目标函数并求解
        object_function = dot(W{t}, phi_x_z{ind}) + sum_cost{ind};
        sol = solvesdp( F{ind}, -object_function, options );

        % 输出得到的各个变量的值
        if sol.problem == 0
            % 这些流程变量的值不是必须的
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

    disp('      更新权向量 w...');
    
    %% 2. 计算 ψ(x,z*,z^)=fai(x,z*)-fai(x,z^) 梯度，求出a和b的值
    % 用 U 来代表这个符号 ψ
    for ind=1:N
        % 梯度 论文中fai(x,z*)-fai(x,z^) 但其 bmrm 代码中用的是 fai(x,z^)-fai(x,z*)
        % 按论文方法似乎不收敛，先用代码方案
        psi_zstar_zhat{ind} = phi_x_z_hat{ind} - phi_x_z_star{ind};
    end
    
    % ======================================================================= %
    % 按照公式计算 At = -1/N* sum( ψ(x,z*,z^) )
    % 按照公式计算 Bt = -1/N* sum( △(z*,z^)+ w'* ψ(x,z*,z^) ) - w'*At
    
    sum_U = zeros(size(W{t}));
    for ind=1:N
        % sum( ψ(x,z*,z^) )
        sum_U = sum_U + psi_zstar_zhat{ind};
        % 保存每一轮中每个样本的损失函数
        sample_loss(t, ind) = delta_zstar_zhat(ind);
    end
    % sum( △(z*,z^) ）
    sum_delta = sum(sample_loss(t,:));
    
    A{t} = 1/N * sum_U;
    B{t} = 1/N *( sum_delta + dot(W{t}, sum_U) ) - dot(W{t}, A{t});
    aver_loss(t) = sum_delta/N;
    ls = aver_loss(t);
    
    fprintf('      平均损失函数△(z*,z^):\t%f\n', aver_loss(t));

    %% 3. 通过求解方程（14）来更新w（包括一个二次规划问题）
	% 将更新后的w保存在 W{t+1}中
    [ kexi W{t+1} obj] = CXSL_Update_W_For_BMRM( A, B, lambda );

    % 计算gap
%     lower_bound = zeros(t,1);
%     for i=1:t
%         lower_bound(i) = A{i}'*W{t+1} + B{i}; % 平均松弛变量的下界
%     end
    % ==================================== %
    % gap 的定义？ 这个地方有问题
    % J代表每一轮的 λΩ(w) + ζ（one slack）（即二次规划的那个目标函数）
    % 直接用损失函数作为gap_cur gap_cur(t+1)表示经过t次循环后的损失函数大小
    gap_cur(t) = aver_loss(t);
    
    % ==================================== %
    % 记录时间
    time(t) = toc;
%     fprintf('      近似间隙 ε:\t%f\n', gap_cur);
    fprintf('      时间花费:\t%1.2f s\n', time(t)); 
   
end

% 循环完成，打印信息
if gap_cur(t) <= gap
    disp('  找到了当前gap下的最优解，算法终止');
    
    % 保存最优分配方案和w
    t_best = t;
    w_best = W{t}; % W{t-1}才是取到最佳 gap 值的那个w，随后更新得到的 W{t} gap可能增大了
    loss_best = gap_cur(t);      
else
    disp('  达到最大循环次数，算法终止');
    t_best = find(aver_loss==min(aver_loss(N+1:end))) %  找到过程中损失最小的那个w作为 w_best
    t_best = t_best(end);
    w_best = W{t_best};
    loss_best = aver_loss(t_best);
end
% 保存最佳w，用于测试其他帧精度
save([ trackpath, '\结构化学习\SSVM_Best_W_New.mat'], 'w_best');

% fprintf('\tw:\t%f\n', w_best);
fprintf('\tt_best:\t%d\n', t_best);
fprintf('\tgap_cur:\t%f\n', loss_best);
fprintf('\ttime consumption:\t%0.2f min\n', sum(time)/60);   
w_for_excel = w_best';

plot(aver_loss, '-*');
% 对得到的收敛曲线进行保存
if 0
    name = 'loss_5_13_initwp_1e-8';
    lossdir = [ trackpath, '\训练结果记录\BMRM\'];
    mkdir(lossdir);
    save([lossdir, name, '.mat'], 'aver_loss','sample_loss','w_best','W');
    saveas(1, [lossdir, name, '.fig']);
end








