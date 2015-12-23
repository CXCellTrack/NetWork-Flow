% ======================================================================= %
%
% FW 标准版 
%
% ======================================================================= %
% 此方法每个样本都保存一个Wi，使得样本更新的时候总是从上一次的Wi=0开始
clear;close all;

[ ~, trackpath ] = getpath( 'training' );
% 载入 CXSL_Test_Linear_all 中计算好的 w 作为初始值（实际发现效果并不好）

load([ trackpath, '\结构化学习\initial_w_New.mat']);
% 注意 w 的顺序不能乱
w = [ wij,bij, wit,bit, wid,bid, wiv,biv, wmj,bmj, wsj,bsj ]'; % 原版b

% 第一个数据集
% w = [-24.05027785	0.965959752	-0.209700235	0.023655542	-0.901444678	0.915527485	-0.723055368	0.78127324	-22.34659216	-3.45491283	-1.682414322	-5.355960441	-2.391659001	2.862181421	-7.382944338	8.382838223	1.94377663	-0.451290137	-1.07738777	-4.844423375	-1.122913059	-0.801496889	3.907101647	-11.61160994	3.710115534	0.998335816	4.252699702	0.790594494	1.207125853	3.799458373	1.390618031	5.18991389	1.129864864	0.673380786	-2.076937813	-1.97433464	-1.980221778	-0.051210814	0.597328997	-3.897482158]';
% 第二个数据集
% w = [-18.21315239	2.048551055	0.090611096	-0.07830978	-0.681768441	0.091705287	-0.284558766	0.113666465	-16.77170209	-2.820207584	-1.606735489	-2.170929556	-2.000511632	2.42450433	-4.406444861	10.34098417	1.758814312	-0.819906672	-2.159095585	-4.969572233	1.29607646	0.202318113	3.379177651	-14.25716614	7.109097539	3.366559674	2.10659084	-2.499899814	-3.9849466	2.501397849	1.351917771	4.699879458	-0.525793863	-0.261628668	-4.945586679	-1.846739083	-4.998199696	-0.003648138	1.737582541	-8.35761154]';
% 第三个数据集
w = [-17.69956928	1.61055833	-0.001787913	0.199592206	-0.287005286	-1.058810184	-0.194715268	0.170821175	-15.86943474	-2.459629861	-1.102503407	-1.233255838	-1.209682446	-2.849135318	-5.932152347	3.635597276	3.112043952	0.905772748	-0.112320949	-0.658275898	4.373948249	0.504445114	-0.211398811	-3.571706443	0	0	0	0	0	0	0	0	0	0	-2.873279295	-0.074622812	-1.438747225	-0.381464188	-2.838608426	-5.896985464]';

if 0
    % initial_w是增广的
    w = zeros(numel(w), 1);
end

% 定义样本个数 N 和 单个样本中的帧数 frame
N = 64;
frame = 2;
s_frame = zeros(N,1);
e_frame = zeros(N,1);
% 目前有gt的帧数，对随机取样有影响
gt_frame = 65;

% 选择取样方式
sample_method = 3;
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
for ii=1:N
    disp(['  ', num2str(s_frame(ii)), '——', num2str(e_frame(ii)), '帧...']);
end

% 初始化： 定义A、B、循环次数上限 iter、间隙阈值 gap
iter = 100;
gap = 0.0010; % 按照 O(1/gap) 的收敛速度，应该在百循环左右完成

%% ============================== 全局变量 =============================== %
gap_cur = zeros(iter,1); % 记录每次得到的gap
gamma = zeros(iter,1); % 步长gamma

W = cell(iter,1); % W 存放综合权值w
W{1} = w; % 全体样本的W需要设定初值w

L = zeros(iter,1); % L 存放综合L
ls = 1; % 样本平均损失函数'

t = 0;
time = zeros(iter,1); % 记录每次循环所用的时间
sample_loss = zeros(iter,N); % 记录每一轮中每个样本的损失函数
aver_loss = zeros(iter,1); % 记录每一轮中样本损失函数均值
% ======================================================================= %
phi_x_z_hat = cell(N,1);
delta_zstar_zhat = zeros(N,1);
U_x_zstar_zhat = cell(N,1);
% ======================================================================== %
% 循环求解部分参数设置
options = sdpsettings('verbose', 0, 'solver', 'gurobi','cachesolvers',1); % cplex设置放到循环外
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

%% 当当前循环次数t小于上限，且gap不符合要求时，进行循环计算，若想增大精度或轮数，修改gap和iter再允许此cell即可
lambda = 1e-2;
usecostall = 0;
linesearch = 0;
if usecostall
    disp('当前选择的损失中包含了虚景！');
    sample_cost = sum_cost_all;
else
    sample_cost = sum_cost;
end

while t < iter % 此处没有采用gap，而是使用loss作为停止判据
    % 记录下每次循环所用的时间
    t = t + 1;
    t0 = clock;
    disp('  ==========================');
    disp(['  开始第 ', num2str(t), ' 轮循环...']);

    %% 1. 计算给定 w 下，每个样本的最佳分配（方程（10））
    
    for ind=1:N
        disp(['      计算样本 ', num2str(ind), '...']); 
        % 临时组建目标函数并求解
        object_function = dot(W{t}, phi_x_z{ind}) + sample_cost{ind};
        sol = solvesdp( F{ind}, -object_function, options );

        % 输出得到的各个变量的值
        if sol.problem == 0      
            phi_x_z_hat{ind} = value(phi_x_z{ind});
            delta_zstar_zhat(ind) = value(sample_cost{ind});
        else
            sol.info
            yalmiperror(sol.problem)
        end
   
    end

    disp('      更新权向量 w...');
    
    %% 2. 计算 ψ(x,z*,z^)=fai(x,z*)-fai(x,z^) 梯度
    % 用 U 来代表这个符号 ψ
    
    sum_U = zeros(size(W{t}));
    for ind=1:N
        % 梯度 论文中fai(x,z*)-fai(x,z^)
        U_x_zstar_zhat{ind} = phi_x_z_star{ind} - phi_x_z_hat{ind};
        % sum( ψ(x,z*,z^) )
        sum_U = sum_U + U_x_zstar_zhat{ind};
        % 保存每一轮中每个样本的损失函数
        sample_loss(t, ind) = delta_zstar_zhat(ind);
    end
    % =================================================================== %
    % sum( △(z*,z^) ）
    sum_delta = sum(sample_loss(t,:));
    aver_loss(t) = sum_delta/N; 
    fprintf('      平均损失函数△(z*,z^):\t%f\n', aver_loss(t));

    %% 3. 通过求解方程（14）来更新w（包括一个二次规划问题）
    % 定义惩罚项 lambda λ，在FW中这个lambda与步长有关，lambda越大，则步长越大
    % 因此当使用svm的w时，由于已经接近标准答案，可使用小lambda微调w即可
	
    Ws = sum_U/(lambda*N);
    ls = aver_loss(t);
    % 下面这几句用于设置初始L，如果注释掉则初始L为0

    % 计算gap，gap的值随样本、lambda都会变化，无法确定下来，因此还是用loss做gap比较合适
    gap_cur(t) = lambda*(W{t}- Ws)'*W{t}- L(t)+ ls;
    % 计算步长gamma
    if linesearch
        gamma(t) = gap_cur(t)/(lambda*norm(W{t}- Ws)^2);
        gamma(t) = max([0, min([gamma(t),1])]);
    else
        gamma(t) = 2*N/(2*N + t-1);
    end
    
    % 更新 w和L，将更新后的w保存在 W{t+1}中
    W{t+1} = (1- gamma(t))*W{t}+ gamma(t)*Ws;
    L(t+1) = (1- gamma(t))*L(t)+ gamma(t)*ls;

    fprintf('      对偶间隙gap:\t%f\n', gap_cur(t));
    % ==================================== %
    % 记录时间
    t1  = clock;
    time(t) = etime(t1, t0);
    fprintf('      时间花费:\t%1.2f s\n', time(t)); 

end

% 循环完成，打印信息
if ls <= gap
    disp('  找到了当前gap下的最优解，算法终止');
    % 保存最优分配方案和w
    t_best = t;
    w_best = Wavg{t}; % W{t}才是取到最佳 gap 值的那个w
    loss_best = ls;      
else
    disp('  达到最大循环次数，算法终止');
    t_best = find(aver_loss==min(aver_loss)) %  找到过程中损失最小的那个w作为 w_best
    t_best = t_best(end);
    w_best = W{t_best};
    loss_best = aver_loss(t_best);
end
% 保存最佳w，用于测试其他帧精度
save([ trackpath, '\结构化学习\SSVM_Best_W_New.mat'], 'w_best');

fprintf('\n\tt_best:\t%d\n', t_best);
fprintf('\tgap_best:\t%f\n', loss_best);
fprintf('\ttime consumption:\t%0.2f min\n', sum(time)/60);   
w_for_excel = w_best';

plot(aver_loss, '-*');
% 对得到的收敛曲线进行保存
if 0
    name = 'loss_64_2_cons35_cost1_init0p_noline';
    lossdir = [ trackpath, '\训练结果记录\FW\'];
    mkdir(lossdir);
    save([lossdir, name, '.mat'], 'time','w','linesearch','use_op_cons','sample_loss','lambda','w_best','W');
    saveas(1, [lossdir, name, '.fig']);
end









