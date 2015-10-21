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
% 载入 CXSL_Test_Linear_all 中计算好的 w 作为初始值（实际发现效果并不好）
if 1
    load([ trackpath, '\结构化学习\initial_w_New.mat']);
    % 注意 w 的顺序不能乱
    w = [ wij, wit, wid, wiv, wmj, wsj ]';
    clear wij wit wid wiv wmj wsj;
else
    % 随机选取w，种子控制 w 可复现
    rng(0); 
    w = rand(42,1);
%     w = zeros(42,1);
end
% 定义样本个数 N 和 单个样本中的帧数 frame
N = 5;
frame = 13;
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

% ============================== 全局变量 =============================== %
% 全部变量最终都被保存下来，局部变量无需保存

% 初始化： 定义A、B、循环次数上限 iter、间隙阈值 gap
iter = 50;
gap = 0.0010; % 按照 O(1/gap) 的收敛速度，应该在百循环左右完成
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
% 计算 fai(x,z)和△(z*,z)，分配好流程变量
for ind=1:N
    disp('  ==========================');
    disp(['  预计算样本',num2str(ind),'的训练数据...']);
    % ----------------------------------------- %
    % 分配各事件流程变量，预先计算好 fai(x,z)和 △(z*,z)
    [ fij{ind} fit{ind} fid{ind} fiv{ind} fmj{ind} fsj{ind} phi_x_z{ind} sum_cost{ind} ] =...
        CXSL_Calculate_Fai_And_Loss( s_frame(ind), e_frame(ind) );
    % 计算约束条件 F，调用 CXSL_Calculate_Constraint_New_Conflict 这个函数
    % 与 BundleMethod_Output_Test 中的同名函数一样
    % ----------------------------------------- %
    % 2015.7.6 使用了新的矛盾约束规则（22矛盾约束）
    [ F{ind} ] = CXSL_Calculate_Constraint_New_Conflict( 'training', true, s_frame(ind), e_frame(ind),...
        fij{ind}, fit{ind}, fid{ind}, fiv{ind}, fmj{ind}, fsj{ind} );
    % ----------------------------------------- %
	% 计算标准答案中的fai(x,z*)
	[ phi_x_z_star{ind} ] = CXSL_Calculate_fai_x_zstar_New( s_frame(ind), e_frame(ind), 'star'); 
    % ----------------------------------------- %
end
toc;

%% 当当前循环次数t小于上限，且gap不符合要求时，进行循环计算，若想增大精度或轮数，修改gap和iter再允许此cell即可

options = sdpsettings('verbose', 0, 'solver', 'cplex', 'saveduals', 0); % cplex设置放到循环外

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
    
    % 定义惩罚项 lambda λ，这个越小，则w越大（收敛速度似乎会加快）
    lambda = 1e-6;
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
    gap_best = gap_cur(t);      
else
    disp('  达到最大循环次数，算法终止');
    t_best = find(aver_loss==min(aver_loss(aver_loss~=0))); %  找到过程中损失最小的那个w作为 w_best
    w_best = W{t_best};
    gap_best = aver_loss(t_best);
end
% 保存最佳w，用于测试其他帧精度
save([ trackpath, '\结构化学习\SSVM_Best_W_New.mat'], 'w_best');

% fprintf('\tw:\t%f\n', w_best);
fprintf('\tt_best:\t%d\n', t_best);
fprintf('\tgap_cur:\t%f\n', gap_best);
fprintf('\ttime consumption:\t%0.2f min\n', sum(time)/60);   
w_for_excel = w_best';

plot(aver_loss, '-*');
% 对得到的收敛曲线进行保存
if 0
    name = 'loss_15_10_y';
    lossdir = [ trackpath, '\训练结果记录\BCFWavg_New\'];
    save([lossdir, name, '.mat'], 'aver_loss','sample_loss','w_best','Wavg');
    saveas(1, [lossdir, name, '.fig']);
end








