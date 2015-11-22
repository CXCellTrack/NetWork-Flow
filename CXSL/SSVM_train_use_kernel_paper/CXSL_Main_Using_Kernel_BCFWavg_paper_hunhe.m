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
load([ trackpath, '\结构化学习\initial_w_New.mat']);
% 注意 w 的顺序不能乱
w = cell(1,6); % 将w也分事件
% w{1} = [wij,bij]'; % 采用svm得到的w
% w{2} = [wit,bit]';
% w{3} = [wid,bid]';
% w{4} = [wiv,biv]';
% w{5} = [wmj,bmj]';
% w{6} = [wsj,bsj]';
w{1} = wij'; % 采用svm得到的w
w{2} = wit';
w{3} = wid';
w{4} = wiv';
w{5} = wmj';
w{6} = wsj';

if 0 % 使用全0w
    w{1} = zeros(size(w{1}));
    w{2} = zeros(size(w{2}));
    w{3} = zeros(size(w{3}));
    w{4} = zeros(size(w{4}));
    w{5} = zeros(size(w{5}));
    w{6} = zeros(size(w{6}));
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
for ii=1:N
    disp(['  ', num2str(s_frame(ii)), '——', num2str(e_frame(ii)), '帧...']);
end

% 定义循环次数上限iter和收敛指标gap
iter = 50;
gap = 0.0010;

%% ------------------ 线性w所用到的变量 ---------------------- %
W = cell(iter,6); % W 存放综合权值w
Wavg = cell(iter,6);
Wi = cell(iter,N); % Wi存放样本权值w
for ev=1:6
    W{1,ev} = w{ev};% 全体样本的W需要设定初值w
    Wavg{1,ev} = w{ev};
end
for i=1:N
    Wi{1,i} = cell(1,6); % Wi 存放每次循环中特定样本更新后的Wi
    for ev=1:6
        Wi{1,i}{ev} = w{ev};% 全体样本的W需要设定初值w
    end 
end

L = zeros(iter,6); % L 存放综合L
Li = cell(iter,N);
for i=1:N
    Li{1,i} = zeros(1,6);
end

% 这个就不要L了
ls = 1; % 样本平均损失函数

%% ------------------ 核函数所用到的变量 ----------------- % 
% 设置6个事件的核函数种类以及参数（可选核函数为linear、poly、rbf、sigmoid）
% 6个事件依次为fij、fit、fid、fiv、fmj、fsj
kernel_type = {'linear','linear','rbf','linear','linear','linear'};
cmd = {'','','-g 0.1','','',''};
% 定义惩罚项 lambda λ（用于控制w的数量级）
lambda = 1e-2*ones(1,6);
islinear = strncmp(kernel_type, 'linear', 6); % 逻辑矩阵，用于指示哪些事件是线性

alpha_i = cell(iter,1); % N个样本，每个样本一个alpha_i
alpha_avg = cell(iter,1);
for tt=1:iter
    alpha_i{tt} = cell(N,6);
    alpha_avg{tt} = cell(N,6);
    for ii=1:N*6
        alpha_i{tt}(ii) = {1};
        alpha_avg{tt}(ii) = {1}; 
    end
end

phi_y_i = cell(N,6);

y_i = cell(iter,1); % 支持向量phi对应的y
for tt=1:iter
    y_i{tt} = cell(N,6);
end


Y_i_star_K_Ytrain = cell(1,6);
Y_i_K_Ytrain = cell(1,6);
phi_x_z_hat = cell(1,6);
n_SV = cell(iter,1);
for tt=1:iter
    n_SV{tt} = zeros(N,6);
end

%% ------------------ 循环中用到的变量 -------------------- %
% ======================================================================= %
% 初始化： 定义A、B、循环次数上限 iter、间隙阈值 gap
gap_cur = zeros(iter,1); % 记录每次得到的gap
gamma = zeros(iter,1); % 步长gamma
t = 0;
time = zeros(iter,1); % 记录每次循环所用的时间
sample_loss = zeros(iter,N); % 记录每一轮中每个样本的损失函数
aver_loss = zeros(iter,1); % 记录每一轮中样本损失函数均值
% ======================================================================= %
% 循环求解部分参数设置
options = sdpsettings('verbose', 0, 'solver', 'gurobi');
rng(0); % 含有随机选择部分，需要设定种子
random = 0; % 变量random作为一个flag，为1时是随机抽样，为0时是按顺序抽样
ind = 0;
% ======================================================================= %
disp('  预计算目标函数和约束条件...');

%% 提前计算好 phi(x,z^) phi(x,z*)和△(z*,z^)，还有约束条件，循环中组装目标函数，再求解
use_op_cons = [3 5];

fij = cell(N,1);
fit = cell(N,1);
fid = cell(N,1);
fiv = cell(N,1);
fmj = cell(N,1);
fsj = cell(N,1);
% 存放循环中求解得到的y
Fij = cell(e_frame(N)-1,1);
Fid = cell(e_frame(N)-1,1);
Fiv = cell(e_frame(N)-1,1);
Fit = cell(e_frame(N)-1,1);
Fsj = cell(e_frame(N),1);
Fmj = cell(e_frame(N),1);

F = cell(N,1);
Fkernel = cell(N,1);
for ii=1:N
    F{ii} = lmi;
    Fkernel{ii} = lmi;
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
        CXSL_Calculate_Event_Fai_And_Loss( s_frame(ii), e_frame(ii) );
    % ----------------------------------------- %
    % 2015.7.6 使用了新的矛盾约束规则（22矛盾约束）
    [ F{ii} ] = CXSL_Calculate_Constraint_New_Conflict( 'training', use_op_cons, s_frame(ii), e_frame(ii),...
        fij{ii}, fit{ii}, fid{ii}, fiv{ii}, fmj{ii}, fsj{ii} );
    % 计算核特有的约束（只是为了减小核运算的规模）
    [ Fkernel{ii} ] = calculate_kernel_constraint( islinear, s_frame(ii), e_frame(ii),...
        fij{ii}, fit{ii}, fid{ii}, fiv{ii}, fmj{ii}, fsj{ii} );
    % 2约束合并
    if 1
        disp('  加入了kernel约束，使答案的全0行强制为0');
        F{ii} = [ F{ii}, Fkernel{ii} ];
    end
    % ----------------------------------------- %
	% 计算标准答案中的phi(x,z*)
	[ phi_x_z_star{ii} ] = CXSL_Calculate_event_fai_x_zstar( s_frame(ii), e_frame(ii), 'star'); 
    % ----------------------------------------- %
end

%% 最耗时的步骤在这里！！！计算所有特征间的核
use_distrabute = 0;
if ~use_distrabute
    kernel_path = [trackpath, '\核训练\kernel_ff_all.mat'];
    if ~exist(kernel_path, 'file')   
        kernel_ff_all = cell(1,6);
        for ev=1:6
            kernel_ff_all{ev} = cell(gt_frame-1);
        end
        s1 = 1; % fe1的开始帧
        e1 = 65; % fe1的结束帧
        s2 = 1; % fe2的开始帧
        e2 = 65; % fe2的结束帧
        for ev=1:6 
            if ~islinear(ev)
                kernel_ff_all{ev} = ssvm_pre_cal_all_kernel_paper(kernel_ff_all{ev}, gt_frame, ev, kernel_type{ev}, cmd{ev}, [s1 e1 s2 e2]);
                disp('保存kernel数据...');tic % 保存这么大的数据也很花时间
                save(kernel_path, 'kernel_ff_all','-v7.3');toc
            end
        end
    else
        disp('载入事先计算好的核函数...');
        tic
        load(kernel_path); % 若存在直接载入即可
        toc
    end

else
    % 使用分散存储
    kernel_ff_all = cell(1,6);
    s1 = 1; % fe1的开始帧
    e1 = 80; % fe1的结束帧
    s2 = 1; % fe2的开始帧
    e2 = 80; % fe2的结束帧
    for ev=1:6 
        ev_kernel_name = [trackpath, '\核训练\dis_kernel_data\row_1_ev',num2str(ev),'.mat'];
        if ~islinear(ev) && ~exist(ev_kernel_name, 'file') 
            ssvm_pre_cal_all_kernel_paper_dis(gt_frame, ev, kernel_type{ev}, cmd{ev}, [s1 e1 s2 e2]);
            for i_file=s1:e1-1
                % 载入分布式存储mat，并组合成大矩阵（慢、卡）
                tmp_path = [trackpath, '\核训练\dis_kernel_data\row_',num2str(i_file),'_ev',num2str(ev),'.mat'];
                disp('载入第',num2str(i_file),'行kernel数据...');
                load(tmp_path);
                kernel_ff_all{ev} = [kernel_ff_all{ev}; rowkernel];
            end
        end
    end
end
    
%% 给y_i、phi_y_i分配初值
for ii=1:N
    for ev=1:6
        if ~islinear(ev) % 只针对非线性核
            % 注意phi_y_i的i是A的索引，而非样本的索引！
            % 初始phi_y设置为标准答案，对应于psi为0，w为全0列
            phi_y_i{ii,ev} = phi_x_z_star{ii}{ev}; 
            % 使用函数来给y_i分配初值
            y_i{1}{ii,ev}{1} = init_assgin_y_i( trackpath, s_frame(ii), e_frame(ii), ev);
        end
    end
end

% global Kernel_ev; % 分配全局变量以加快速度(观测速度无明显提升，还是不用了)

%% 当当前循环次数t小于上限，且gap不符合要求时，进行循环计算，若想增大精度或轮数，修改gap和iter再运行此cell即可
usecostall = 0;
linesearch = 1;
if usecostall
    disp('当前选择的损失中包含了虚景！');
    sample_cost = sum_cost_all;
else
    sample_cost = sum_cost;
end

while t < iter %&& ls*N >= gap) || t <= N % 迭代次数必须大于样本数（即每个样本都必须用到）
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
    disp(['      计算样本 ', num2str(ind), '...']); 
    % 计算总样本的alpha和phi_y
    for ii=1:N
        for ev=1:6
            % 计算各个样本的支持向量个数（可供查看）
            n_SV{t}(ii,ev) = numel(alpha_i{t}{ii,ev});
        end
    end

    %% 2. 临时组建目标函数并求解
    % 目标函数表达式：
    K_OBJ = 0;
    W_phi_all = 0;
    tic
    for ev=1:6
        if ~islinear(ev) % 只针对非线性核
            %% ------------- 非线性核目标函数 ----------------- %
            % 使用非线性核时这个内积用其他核函数代替
            % 并且按不同事件可以使用不同核
%             profile on
            Kernel_ev = kernel_ff_all{ev};
            y_ev = y_i{t}(:,ev); % 取到该事件的y_i
            alpha_ev = alpha_avg{t}(:,ev); % 取到该事件的alpha_i
            ind_train = ind;
            % 1、计算<phi_i*,phi>
            % 代表Y_i'*K*Y 
            % 2、计算<phi_i,phi>
            % 代表Y_i'*K*Y
            [ Y_i_star_K_Ytrain{ev} Y_i_K_Ytrain{ev} ] = ssvm_kernel_train_paper(y_ev, alpha_ev, Kernel_ev, ev, N, ind_train, s_frame, e_frame,...
                                            fij, fit, fid, fiv, fmj, fsj); 
            % 3、合成目标函数
            K_OBJ = K_OBJ + (Y_i_star_K_Ytrain{ev} - Y_i_K_Ytrain{ev})/(lambda(ev)*N); % 目标函数表达式需要将所有事件加起来
%             profile viewer
            
        else
            %% --------------- 线性核目标函数 ---------------- %
            % 对于线性来说，只需导出<w,phi>即可
            W_phi_all = W_phi_all + dot(Wavg{t,ev}, phi_x_z{ind}{ev});
        
        end
    end
    toc;
    object_function = K_OBJ + W_phi_all + sample_cost{ind};
    disp('    目标函数组建完毕');
    sol = solvesdp( F{ind}, -object_function, options );

    %% 3. 输出得到的各个变量的值（即y_i，可认为是支持向量）
    if sol.problem == 0      
        for ev=1:6
            phi_x_z_hat{ev} = value(phi_x_z{ind}{ev});
        end
        delta_zstar_zhat = value(sample_cost{ind});
        
        % 需要把phi_x_z_hat对应的那个Y也求处来
        for ff = s_frame(ind):e_frame(ind)-1
            Fij{ff} = round(value(fij{ind}{ff})) ;
            Fid{ff} = round(value(fid{ind}{ff})) ;
            Fiv{ff} = round(value(fiv{ind}{ff})) ;
            Fit{ff} = round(value(fit{ind}{ff})) ;
        end
        for ff = s_frame(ind)+1:e_frame(ind)
            Fsj{ff} = round(value(fsj{ind}{ff})) ;
            Fmj{ff} = round(value(fmj{ind}{ff})) ;
        end
        Fsj(s_frame(ind)) = []; % 去除第一位，保持和其他事件一样长度（但实际上往前推了一帧）
        Fmj(s_frame(ind)) = []; % 去除第一位，保持和其他事件一样长度（但实际上往前推了一帧）
        
    else
        sol.info
        yalmiperror(sol.problem)
    end
   
    disp('      更新对偶变量 α...');
    
    %% 4. 判断s对应的phi之前是否出现过，并更新alpha_i和y_i
    % 将得到的phi_x_z_hat与之前所有的phi进行比较
    % 实际上，如果2个y不相同，但得到的phi有可能相同，因此通过比较phi来确定y是否一样不严谨
    % 正确的方法应该通过比较y来进行，但只要phi相同，alpha对phi组合之后的结果还是一样，因此无本质区别
    %
    % 将其他的样本从上一轮带过来，再更新ind样本
    alpha_i{t+1} = alpha_i{t};
    alpha_avg{t+1} = alpha_avg{t}; % 这2局如果写在事件循环里面会出错！
    %
    
    for ev=1:6 % 每个事件分别更新
        %% --------------- 线性核更新方法 ---------------- %
        if islinear(ev)
            PSI_zstar_zhat = phi_x_z_star{ind}{ev} - phi_x_z_hat{ev};
            Ws = PSI_zstar_zhat/(lambda(ev)*N);
            ls = delta_zstar_zhat/N;  
            
            if linesearch
                % line-search寻找最佳gamma
                tmp = lambda(ev)*(Wi{t,ind}{ev}- Ws)'*W{t,ev}- Li{t,ind}(ev)+ ls;
                gamma(t) = tmp/(lambda(ev)*norm(Wi{t,ind}{ev}- Ws)^2);
                gamma(t) = max([0, min([gamma(t),1])]); % clip to 0-1
            else
                % 普通方法计算步长gamma
                gamma(t) = 2*N/(2*N + t-1);
            end
            % 更新 wi和Li，将更新后的w保存在 W{t+1,ind}中
            Wi{t+1,ind}{ev} = (1- gamma(t))*Wi{t,ind}{ev} + gamma(t)*Ws;
            Li{t+1,ind}(ev) = (1- gamma(t))*Li{t,ind}(ev) + gamma(t)*ls;
            % 对于此轮没轮到的样本，将其Wi和Li带入下一轮中
            sample_not_used = mysetdiff(1:N, ind);
            for jj=sample_not_used
                % jj 没被用到的样本编号
                Wi{t+1,jj}{ev} = Wi{t,jj}{ev}; % 直接将其Wi带入下一轮
                Li{t+1,jj}(ev) = Li{t,jj}(ev);
            end

            % 更新 w和L，将更新后的w保存在 W{t+1,N+1}中
            W{t+1,ev} = W{t,ev} + Wi{t+1,ind}{ev} - Wi{t,ind}{ev};
            
            % 更新Wavg
            Wavg{t+1,ev} = (t-1)/(t+1)*Wavg{t,ev} + 2/(t+1)*W{t+1,ev};
        end
        
        %% --------------- 非线性核更新方法 ---------------- %
        if ~islinear(ev)
            flag_equal = 0;
            n_phi = size(phi_y_i{ind,ev},2); % 当前样本支持向量的个数
            for m=1:n_phi
                if isequal(phi_y_i{ind,ev}(:,m), phi_x_z_hat{ev}) 
                    flag_equal = true;
                    break;
                end
            end
            % ------------------------------------------- %
            if flag_equal
                % 如果有重复，则新一轮的phi不变
                phi_y_i{ind,ev} = phi_y_i{ind,ev};
                s = m; % 找出s的位置
                n_phi_new = n_phi; 
                y_i{t+1} = y_i{t}; % y_i不变
            else
                phi_y_i{ind,ev} = [ phi_y_i{ind,ev}, phi_x_z_hat{ev}]; % 否则将此轮得到的phi加入新一轮中
                s = size(phi_y_i{ind,ev},2); % s为新出现的
                n_phi_new = n_phi + 1; % 更新后的支持向量数目
                % 同时其对应的y^也要加入到y_i中
                y_i{t+1} = y_i{t}; % 将其他的样本从上一轮带过来，再更新ind样本
                y_i{t+1}{ind,ev}{n_phi_new} = assign_y_hat_into_y_i(ev, ind, s_frame, e_frame,Fij,Fit,Fid,Fiv,Fmj,Fsj);
                
            end

            alpha_vector = zeros(n_phi_new,1);
            alpha_vector(1:n_phi) = alpha_i{t}{ind,ev}; % 旧aplha向量
            s_vector = zeros(n_phi_new,1);
            s_vector(s) = 1; % 新aplha向量
            % 更新alpha_i
            gamma_alpha = 2*N/(2*N + t-1);
            alpha_i{t+1}{ind,ev} = (1-gamma_alpha)*alpha_vector + s_vector*gamma_alpha;

            % 更新alpha_avg（向量长度不等则增加一位，变为等长）
            if numel(alpha_avg{t}{ind,ev})~=numel(alpha_i{t+1}{ind,ev})
                alpha_avg{t}{ind,ev} = [alpha_avg{t}{ind,ev}; 0];
            end
            alpha_avg{t+1}{ind,ev} = (t-1)/(t+1)*alpha_avg{t}{ind,ev} + 2/(t+1)*alpha_i{t+1}{ind,ev};
            
        end
        
    end

    %% 5. 统计损失、时间花费等数据
    % =================================================================== %
    sample_loss(t, ind) = delta_zstar_zhat;
    aver_loss(t) = delta_zstar_zhat; 
    fprintf('      当前样本损失函数△(z*,z^):\t%f\n', aver_loss(t));
    % 记录时间
    time(t) = etime(clock, tstart);
    fprintf('      时间花费:\t%1.2f s\n', time(t)); 
    
end

% 循环完成，打印信息
if ls*N <= gap
    disp('  找到了当前gap下的最优解，算法终止');
    % 保存最优分配方案和w
    t_best = t;
else
    disp('  达到最大循环次数，算法终止');
    t_best = find(aver_loss==min(aver_loss(N+1:end))) %  找到过程中损失最小的那个w作为 w_best
    t_best = t_best(t_best>N);
    t_best = t_best(end);
end   

sample_id = mod(t_best,N); % 根据tbest得到此时的样本编号
sample_id(sample_id==0) = N;

%% 得到最终的y^，y*和alpha_i，在测试中使用这2个即可
w_best = Wavg(t_best,:);
yhat_best = y_i{t_best};

load([ trackpath, '\结构化学习\Feature_Plus_New.mat']); % 载入特征
load([ trackpath, '\GT\GT_after_hand_tune\GT_Flow_Variables_New.mat']); % 载入标准答案
ystar_best = cell(N,6);
feature_best = cell(N,6);
for ii=1:N % 求标准答案y*
    ystar_best{ii,1} = Fij(s_frame(ii):e_frame(ii)-1);
    ystar_best{ii,2} = Fit(s_frame(ii):e_frame(ii)-1);
    ystar_best{ii,3} = Fid(s_frame(ii):e_frame(ii)-1);
    ystar_best{ii,4} = Fiv(s_frame(ii):e_frame(ii)-1);
    ystar_best{ii,5} = Fmj(s_frame(ii)+1:e_frame(ii));
    ystar_best{ii,6} = Fsj(s_frame(ii)+1:e_frame(ii));
    
    % 找出对应的特征
    feature_best{ii,1} = feature_fij_p(s_frame(ii):e_frame(ii)-1);
    feature_best{ii,2} = feature_fit_p(s_frame(ii):e_frame(ii)-1);
    feature_best{ii,3} = feature_fid_p(s_frame(ii):e_frame(ii)-1);
    feature_best{ii,4} = feature_fiv_p(s_frame(ii):e_frame(ii)-1);
    feature_best{ii,5} = feature_fmj_p(s_frame(ii)+1:e_frame(ii));
    feature_best{ii,6} = feature_fsj_p(s_frame(ii)+1:e_frame(ii));
end
% 减少一下数据的规模
for ev=1:6 
    if islinear(ev)
        ystar_best(:,ev) = {[]};
        feature_best(:,ev) = {[]};
    end
end

alpha_best = alpha_avg{t_best};
n_SV_best = n_SV{t_best};

% 最终使用 y*-y^ 进行预测
delta_y_best = cell(N,6);
for nn=1:N*6
    delta_y_best{nn} = {};
end
for ev=1:6 % 事件循环，只招核
    if islinear(ev)
        continue;
    end
    for ii=1:N % 样本循环
        for nsv=1:n_SV_best(ii,ev) % 每个样本中的支持向量数目
            for tt=1:numel(ystar_best{ii,ev}) % 每个样本帧数
                % 即y*-y^（测试时这2个都是定值，可以先减）
                delta_y_best{ii,ev}{nsv,tt} = ystar_best{ii,ev}{tt} - yhat_best{ii,ev}{nsv}{tt};
                % 第ii个样本第ev个事件中，第nsv个alpha(ii)对应的第tt帧
            end
        end
    end
end

loss_best = aver_loss(t_best);  

%% 保存最佳w，用于测试其他帧精度

fprintf('\n\tt_best:\t%d\n', t_best);
fprintf('\tgap_best:\t%f\n', loss_best);
fprintf('\ttime consumption:\t%0.2f min\n', sum(time)/60);   

plot(aver_loss, '-*');
% 对得到的收敛曲线进行保存
if 0
    name = 'loss_5_13_initwp_line';
    lossdir = [ trackpath, '\训练结果记录\核记录\'];
    mkdir(lossdir);
    save([lossdir, name, '.mat'], 'time','sample_loss','w_best','linesearch',...
        'delta_y_best','feature_best','alpha_best','kernel_type','cmd','lambda','islinear');
    saveas(1, [lossdir, name, '.fig']);
end







