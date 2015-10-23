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
clear;close all;

[ ~, trackpath ] = getpath( 'training' );
load([ trackpath, '\结构化学习\Feature_New.mat']); % 载入特征
load([ trackpath, '\GT\GT_after_hand_tune\GT_Flow_Variables_New.mat']); % 载入标准答案

% 定义样本个数 N 和 单个样本中的帧数 frame
N = 3;
frame = 5;
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

%% 定义循环次数上限iter和收敛指标gap
iter = 50;
gap = 0.0010;

%% ------------------ 核函数所用到的变量 ----------------- % 
% 设置6个事件的核函数种类以及参数（可选核函数为linear、poly、rbf、sigmoid）
% 6个事件依次为fij、fit、fid、fiv、fmj、fsj
kernel_type = {'linear','linear','linear','linear','linear','linear'};
cmd = {'-d 2 -g 1','-d 2 -g 1','-d 2 -g 1','-d 2 -g 1','-d 2 -g 1','-d 2 -g 1'};

alpha_i = cell(N,6); % N个样本，每个样本一个alpha_i
for ii=1:N
    for ev=1:6
        alpha_i{ii,ev} = 1;
    end
end

phi_y_i = cell(N,6);
y_i = cell(N,6); % 支持向量phi对应的y

K_phi = cell(1,6);
phi_x_z_hat = cell(1,6);
n_SV = zeros(N,6);

%% ------------------ 循环中用到的变量 -------------------- %
% ======================================================================= %
% 初始化： 定义A、B、循环次数上限 iter、间隙阈值 gap
gap_cur = zeros(iter,1); % 记录每次得到的gap
gamma = zeros(iter,1); % 步长gamma
ls = 1;
t = 0;
time = zeros(iter,1); % 记录每次循环所用的时间
sample_loss = zeros(iter,N); % 记录每一轮中每个样本的损失函数
aver_loss = zeros(iter,1); % 记录每一轮中样本损失函数均值
% ======================================================================= %
% 循环求解部分参数设置
% options = sdpsettings('verbose', 0, 'solver', 'gurobi');
options = sdpsettings('verbose', 0, 'solver', 'cplex', 'saveduals', 0); % cplex设置放到循环外
rng(0); % 含有随机选择部分，需要设定种子
random = 0; % 变量random作为一个flag，为1时是随机抽样，为0时是按顺序抽样
ind = 0;
% ======================================================================= %
disp('  预计算目标函数和约束条件...');

%% 提前计算好 phi(x,z^) phi(x,z*)和△(z*,z^)，还有约束条件，循环中组装目标函数，再求解
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
phi_x_z_star = cell(N,1);

tic;
% 计算 phi(x,z)和△(z*,z)，分配好流程变量
for ii=1:N
    disp('  ==========================');
    disp(['  预计算样本',num2str(ii),'的训练数据...']);
    % ----------------------------------------- %
    % 分配各事件流程变量，预先计算好 phi(x,z)和 △(z*,z)
    [ fij{ii} fit{ii} fid{ii} fiv{ii} fmj{ii} fsj{ii} phi_x_z{ii} sum_cost{ii} ] =...
        CXSL_Calculate_Event_Fai_And_Loss( s_frame(ii), e_frame(ii) );
    % 计算约束条件 F，调用 CXSL_Calculate_Constraint_New_Conflict 这个函数
    % 与 BundleMethod_Output_Test 中的同名函数一样
    % ----------------------------------------- %
    % 2015.7.6 使用了新的矛盾约束规则（22矛盾约束）
    [ F{ii} ] = CXSL_Calculate_Constraint_New_Conflict( 'training', true, s_frame(ii), e_frame(ii),...
        fij{ii}, fit{ii}, fid{ii}, fiv{ii}, fmj{ii}, fsj{ii} );
    % ----------------------------------------- %
	% 计算标准答案中的phi(x,z*)
	[ phi_x_z_star{ii} ] = CXSL_Calculate_event_fai_x_zstar( s_frame(ii), e_frame(ii), 'star'); 
    % ----------------------------------------- %
end
toc;

%% 当当前循环次数t小于上限，且gap不符合要求时，进行循环计算，若想增大精度或轮数，修改gap和iter再运行此cell即可
for ii=1:N
    for ev=1:6
        % 注意phi_y_i的i是A的索引，而非样本的索引！
        % 初始phi_y设置为标准答案，对应于psi为0，w为全0列
        phi_y_i{ii,ev} = phi_x_z_star{ii}{ev}; 
        % 使用函数来给y_i分配处值
        y_i{ii,ev}{1} = init_assgin_y_i( trackpath, s_frame(ii), e_frame(ii), ev);
    end
end

while (t < iter && ls >= gap) || t <= N % 迭代次数必须大于样本数（即每个样本都必须用到）
    t = t + 1;
    
    % 记录下每次循环所用的时间
    tstart = clock;
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
    
    %% 1. 计算给定 aplha 下，当前选定样本的最佳分配（方程（10））

    % 定义惩罚项 lambda λ
    lambda = 1e-2;
    disp(['      计算样本 ', num2str(ind), '...']); 
    % 计算总样本的alpha和phi_y
    for ii=1:N
        for ev=1:6
            % 计算各个样本的支持向量个数（可供查看）
            n_SV(ii,ev) = numel(alpha_i{ii,ev});
        end
    end

    %% 2. 临时组建目标函数并求解
    % 目标函数表达式：
    % H(y） = sigma[ alpha(i)*K(psi(i), phi(y)) ] + L(y)*lambda*N
    % 将y拆出来后变为
    % H(y） = sigma[ alpha(i)*y(i)*K(psi(i), feature(i)) ] + L(y)*lambda*N
    % ------------------------------------- %
    Predict_OBJ = 0;
    for ev=1:6
        % 使用非线性核时这个内积用其他核函数代替
        % 并且按不同事件可以使用不同核
        
        % 使用info来记录phi2的相关信息（计算kernel是用的feature而不是phi2）
        info.ind = ind;
        info.N = N;
        info.ev = ev;
        info.s_frame = s_frame(ind);
        info.e_frame = e_frame(ind);
        % 选择核（目前有2种feature方案，加1的增广和不加的原始）
        % 对于某样本的某事件，进行kernel运算
        %% 1、计算<phi_i*,phi>
        K_phi_i_stat{ev} = ssvm_kernel_train_paper(n_SV, alpha_i, info, kernel_type{ev}, cmd{ev},...
                    fij, fit, fid, fiv, fmj, fsj, s_frame, e_frame,...
                    feature_fij, feature_fit, feature_fid, feature_fiv, feature_fmj, feature_fsj);   
                
        
        %% 2、计算<phi_i,phi>
        y_kernel = y_i(:,ev); % 取到该事件的y_i
        K_phi_i{ev} = ssvm_kernel_train_paper(n_SV, alpha_i, info, kernel_type{ev}, cmd{ev},...
                    fij, fit, fid, fiv, fmj, fsj, s_frame, e_frame,...
                    feature_fij, feature_fit, feature_fid, feature_fiv, feature_fmj, feature_fsj);
        
        Predict_OBJ = Predict_OBJ + K_phi_i_stat{ev} - K_phi_i{ev}; % 目标函数表达式需要将所有事件加起来
    end
    object_function = Predict_OBJ + sum_cost{ind}*lambda*N;
    sol = solvesdp( F{ind}, -object_function, options );

    % 输出得到的各个变量的值
    if sol.problem == 0      
        for ev=1:6
            phi_x_z_hat{ev} = value(phi_x_z{ind}{ev});
        end
        delta_zstar_zhat = value(sum_cost{ind});
    else
        sol.info
        yalmiperror(sol.problem)
    end
   
    disp('      更新对偶变量 α...');
    
    %% 3. 判断s对应的phi之前是否出现过，并更新alpha
    % 将得到的phi_x_z_hat与之前所有的phi进行比较
    % 实际上，如果2个y不相同，但得到的phi有可能相同，因此通过比较phi来确定y是否一样不严谨
    % 正确的方法应该通过比较y来进行，但只要phi相同，alpha对phi组合之后的结果还是一样，因此无本质区别
    for ev=1:6 % 每个事件分别更新
        flag_equal = 0;
        n_phi = size(phi_y_i{ind}{ev},2); % 当前样本支持向量的个数
        for m=1:n_phi
            if isequal(phi_y_i{ind}{ev}(:,m), phi_x_z_hat{ev}) 
                flag_equal = true;
                break;
            end
        end
        % ------------------------------------------- %
        if flag_equal
            % 如果有重复，则新一轮的phi不变
            phi_y_i{ind}{ev} = phi_y_i{ind}{ev};
            s = m; % 找出s的位置
            n_phi_new = n_phi; 
        else
            phi_y_i{ind}{ev} = [ phi_y_i{ind}{ev}, phi_x_z_hat{ev}]; % 否则将此轮得到的phi加入新一轮中
            s = size(phi_y_i{ind}{ev},2); % s为新出现的
            n_phi_new = n_phi + 1; % 更新后的支持向量数目
        end

        alpha_vector = zeros(n_phi_new,1);
        alpha_vector(1:n_phi) = alpha_i{ind}{ev};
        s_vector = zeros(n_phi_new,1);
        s_vector(s) = 1; 
        % 更新alpha_i
        gamma = 2*N/(2*N + t-1);
        alpha_i{ind}{ev} = (1-gamma)*alpha_vector + s_vector*gamma;
    end

    %% 统计损失、时间花费等数据
    % =================================================================== %
    sample_loss(t, ind) = delta_zstar_zhat;
    ls = sum(sample_loss(t,:));
    aver_loss(t) = ls; 
    fprintf('      当前样本损失函数△(z*,z^):\t%f\n', aver_loss(t));
    % 记录时间
    time(t) = etime(clock, tstart);
%     fprintf('      近似间隙 ε:\t%f\n', gap_cur);
    fprintf('      时间花费:\t%1.2f s\n', time(t)); 
    
end

% 循环完成，打印信息
if ls*N <= gap
    disp('  找到了当前gap下的最优解，算法终止');
    % 保存最优分配方案和w
    t_best = t;
else
    disp('  达到最大循环次数，算法终止');
    t_best = find(aver_loss==min(aver_loss(aver_loss~=0))); %  找到过程中损失最小的那个w作为 w_best   
    t_best
    t_best = t_best(end);
end   

sample_id = mod(t_best,N); % 根据tbest得到此时的样本编号
sample_id(sample_id==0) = N;

% 得到最终的A和alpha，在测试中使用这2个即可
A_best = phi_y_all{t_best};
alpha_best = alpha_all{t_best};
loss_best = aver_loss(t_best);  

% 保存最佳w，用于测试其他帧精度
save([ trackpath, '\结构化学习\SSVM_Best_W_New.mat'], 'A_best','alpha_best','kernel_type','cmd');

fprintf('\n\tt_best:\t%d\n', t_best);
fprintf('\tgap_best:\t%f\n', loss_best);
fprintf('\ttime consumption:\t%0.2f min\n', sum(time)/60);   

plot(aver_loss, '-*');
% 对得到的收敛曲线进行保存
if 0
    name = 'loss_15_10_y';
    lossdir = [ trackpath, '\训练结果记录\BCFW_kernel\'];
    save([lossdir, name, '.mat'], 'sample_loss','A_best','alpha_best','kernel_type','cmd');
    saveas(1, [lossdir, name, '.fig']);
end







