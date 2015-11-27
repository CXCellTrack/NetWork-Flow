% ======================================================================= %
%
% BCFWavg_my 与原始的BCFWavg略有不同，原理不太对，但有时效果比较好！
%
% ======================================================================= %
clear;close all;

[ ~, trackpath ] = getpath( 'training' );
% 载入 CXSL_Test_Linear_all 中计算好的 w 作为初始值（实际发现效果并不好）
load([ trackpath, '\结构化学习\initial_w_New.mat']);
% 注意 w 的顺序不能乱
if 0
    w = [ wij, wit, wid, wiv, wmj, wsj ]';
else
%     w = [ wij,bij, wit,bit, wid,bid, wiv,biv, wmj,bmj, wsj,bsj ]'; % 原版b
    % 第二个数据集
%     w = [-18.21315239	2.048551055	0.090611096	-0.07830978	-0.681768441	0.091705287	-0.284558766	0.113666465	-16.77170209	-2.820207584	-1.606735489	-2.170929556	-2.000511632	2.42450433	-4.406444861	10.34098417	1.758814312	-0.819906672	-2.159095585	-4.969572233	1.29607646	0.202318113	3.379177651	-14.25716614	7.109097539	3.366559674	2.10659084	-2.499899814	-3.9849466	2.501397849	1.351917771	4.699879458	-0.525793863	-0.261628668	-4.945586679	-1.846739083	-4.998199696	-0.003648138	1.737582541	-8.35761154]';
    % 第三个数据集
    w = [-17.69956928	1.61055833	-0.001787913	0.199592206	-0.287005286	-1.058810184	-0.194715268	0.170821175	-15.86943474	-2.459629861	-1.102503407	-1.233255838	-1.209682446	-2.849135318	-5.932152347	3.635597276	3.112043952	0.905772748	-0.112320949	-0.658275898	4.373948249	0.504445114	-0.211398811	-3.571706443	0	0	0	0	0	0	0	0	0	0	-2.873279295	-0.074622812	-1.438747225	-0.381464188	-2.838608426	-5.896985464]';
end
clear wij wit wid wiv wmj wsj;
% 也可以选择随机的w或全0的w
if 0
    % initial_w是增广的
    w = zeros(numel(w), 1);
end

% 定义样本个数 N 和 单个样本中的帧数 frame
N = 8;
frame = 8;
s_frame = zeros(N,1);
e_frame = zeros(N,1);
% 目前有gt的帧数，对随机取样有影响
gt_frame = 65;

% 选择取样方式
% 1为接龙取样，2为滑窗取样，3为随机取样
sample_method = 1;
switch sample_method
    case 1 % 取样方法1：接龙取样
        % 按 1-5, 6-10, 11-15, 16-20 这样的方法取样本
        for ind=1:N
            s_frame(ind) = (ind - 1)*frame + 1;
            e_frame(ind) = s_frame(ind) + frame - 1;
        end

    case 2 % 取样方法2：滑窗取样
        % 按 1-5，2-6，3-7 这样的方法取样本
        for ind=1:N
            s_frame(ind) = ind;
            e_frame(ind) = s_frame(ind) + frame - 1;
        end
       
    case 3 % 取样方法3：随机取样
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
iter = 200;
gap = 0.0010; % 按照 O(1/gap) 的收敛速度，应该在百循环左右完成

%% ============================== 全局变量 =============================== %
gap_cur = zeros(iter,1); % 记录每次得到的gap
gamma = zeros(iter,1); % 步长gamma

W = cell(iter,1); % W 存放综合权值w
Wi = w; % Wi 存放每次循环中特定样本更新后的Wi
Wi_old = 0; % Wi_old 存放每次循环中特定样本更新前的Wi
Wavg = cell(iter,1); % Wavg
Ws = 0;

L = zeros(iter,1); % L 存放综合L
Li = 0; % Li 存放每次选中的样本更新后的L
Li_old = 0; % Li_old 存放每次选中的样本更新前的L
ls = 1; % 样本平均损失函数

W{1} = w; % 全体样本的W需要设定初值w
Wavg{1} = w; % Wavg设置初值也为w

t = 0;
time = zeros(iter,1); % 记录每次循环所用的时间
sample_loss = zeros(iter,N); % 记录每一轮中每个样本的损失函数
aver_loss = zeros(iter,1); % 记录每一轮中样本损失函数均值
% ======================================================================= %
% 循环求解部分参数设置
options = sdpsettings('verbose', 0, 'solver', 'gurobi'); % cplex设置放到循环外
rng(0); % 含有随机选择部分，需要设定种子
random = 0; % 变量random作为一个flag，为1时是随机抽样，为0时是按顺序抽样
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
    [ fij{ii} fit{ii} fid{ii} fiv{ii} fmj{ii} fsj{ii} phi_x_z{ii} sum_cost{ii} ] =...
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
lambda = 1e-2;
linesearch = 1;
usecostall = 0;
if usecostall
    sum_cost = sum_cost_all;
end

while t < iter %&& ls*N >= gap || t <= N % 迭代次数必须大于样本数（即每个样本都必须用到）
    t = t + 1;
    
    % 记录下每次循环所用的时间
    tic;
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
    object_function = dot(Wavg{t}, phi_x_z{ind}) + sum_cost{ind};
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
    % 梯度 论文中phi(x,z*)-phi(x,z^)
    U_x_zstar_zhat = phi_x_z_star{ind} - phi_x_z_hat;
    % sum( ψ(x,z*,z^) )
    sum_U = U_x_zstar_zhat;
    % 保存每一轮中每个样本的损失函数
    sample_loss(t, ind) = delta_zstar_zhat;
    % =================================================================== %
    % sum( △(z*,z^) ）
    sum_delta = sum(sample_loss(t,:));
    aver_loss(t) = sum_delta/1; 
    fprintf('      当前样本损失函数△(z*,z^):\t%f\n', aver_loss(t));

    %% 3. 求解最优步长来更新Wavg
    Ws = sum_U/(lambda*N);
    ls = sum_delta/N;
    
    % 计算gap，gap的值随样本、lambda都会变化，无法确定下来，因此还是用loss做gap比较合适
    % 计算步长gamma
    if linesearch
        gap_cur(t) = lambda*(Wi- Ws)'*W{t}- Li+ ls;
        gamma(t) = gap_cur(t)/(lambda*norm(Wi- Ws)^2); % line search
        gamma(t) = max([0, min([gamma(t),1])]);
    else
        gamma(t) = 2*N/(2*N + t-1); % 普通方法
    end
    % 更新 wi和Li，将更新后的w保存在 W{t+1,ind}中
    Wi_old = Wi;
    Li_old = Li;
    Wi = (1- gamma(t))*Wi_old+ gamma(t)*Ws;
    Li = (1- gamma(t))*Li_old+ gamma(t)*ls;
    
    % 更新 w和L，将更新后的w保存在 W{t+1,N+1}中
    W{t+1} = W{t} + Wi - Wi_old;
    L(t+1) = L(t) + Li - Li_old; % 计算L似乎没什么用？
    
    % 更新Wavg
    Wavg{t+1} = (t-1)/(t+1)*Wavg{t} + 2/(t+1)*W{t+1};

    fprintf('      对偶间隙gap:\t%f\n', gap_cur(t));
    % ==================================== %
    % 记录时间
    time(t) = toc;
%     fprintf('      近似间隙 ε:\t%f\n', gap_cur);
    fprintf('      时间花费:\t%1.2f s\n', time(t)); 

end

% 循环完成，打印信息
if ls*N <= gap
    disp('  找到了当前gap下的最优解，算法终止');
    % 保存最优分配方案和w
    t_best = t;
    w_best = Wavg{t_best}; % W{t}才是取到最佳 gap 值的那个w
    gap_best = ls*N;      
else
    disp('  达到最大循环次数，算法终止');
    t_best = find(aver_loss==min(aver_loss(N+1:end))) %  找到过程中损失最小的那个w作为 w_best
    t_best = t_best(t_best>N);
    t_best = t_best(1);
    w_best = Wavg{t_best};
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
% 对得到的收敛曲线进行保存
if 0
    name = 'loss_8_8_cons35_cost1_initwp_line_b';
    lossdir = [ trackpath, '\训练结果记录\BCFWavg_my\initwp\8_8\'];
    mkdir(lossdir);
    save([lossdir, name, '.mat'], 'w','linesearch','use_op_cons','sample_loss','lambda','w_best','Wavg');
    saveas(1, [lossdir, name, '.fig']);
end







