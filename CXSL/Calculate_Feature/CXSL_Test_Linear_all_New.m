%
% ======================================================= %
%
% 使用libsvm测试特征是否线性可分
% 这个函数运用全体特征进行测试，与 CXSL_Test_Linear 有所区别
% 将计算出来较好分类的各事件 w(增广) 保存下来，作为结构化学习的初始值
% 
% ======================================================= %

clear;close all;

dataset = 'training'; % 这个脚本一般只在训练集中使用
[ segpath trackpath ] = getpath( dataset );

load([ trackpath, '\Pair\Pre_data_New.mat'], 'Ellipse','candidate_k_last','candidate_k_next','n','candidate_fij');
% 载入标准答案GT
load([ trackpath, '\GT\GT_after_hand_tune\GT_Flow_Variables_New.mat']);
% 载入特征（不包含增广）
load([ trackpath, '\结构化学习\Feature_New.mat']);
frame = numel(Fmj);    

%% 提取分裂divide样本

positive_sample = [];
negative_sample = [];
num_p_s = 1;
num_n_s = 1;

for t=1:frame-1
    for j=1:n(t)
        if sum( Fid{t}(j,:) )~=0 % 出现分裂则记下正样本
            mm_p = find(Fid{t}(j,:)==1);
            % 样本储存格式
            positive_sample(num_p_s,:) = [ t, j, mm_p, candidate_k_next{t}{j,mm_p}];
            num_p_s = num_p_s + 1;
            
            % 同时选取其他几个中的一个作为负样本
            mm_ct = find(Fid{t}(j,:)==0);

            mm_n = mm_ct( randi(numel(mm_ct)) );  % 随机选取一个负样本
            negative_sample(num_n_s,:) = [ t, j, mm_n, candidate_k_next{t}{j,mm_n}];
            num_n_s = num_n_s + 1;
            
        else
            % =======================================
            % 设置一个随机过滤器，来减少负样本的采样数量
            if randi(10)~=1
                continue; 
            end % 只有1/10的概率会选择一个负样本
            % =======================================

            mm_n = randi(6);       % 没出现分裂则随机选取一个负样本
            negative_sample(num_n_s,:) = [ t, j, mm_n, candidate_k_next{t}{j,mm_n}];
            num_n_s = num_n_s + 1;
        end
    end
end
num_p_s = num_p_s -1; % 正样本数目
num_n_s = num_n_s -1; % 负样本数目

% 输入到SVM分类器中验证是否线性可分
% 选择 'margin' 模式的时候2/3是失效的，训练样本数目就是边界点数
model = CXSL_SVM( positive_sample, negative_sample, feature_fid, 2/3, 3, 'rand');% 倒数第二位表示测试中负样本数目是正样本的多少倍

coef = repmat( model.sv_coef, 1, size(model.SVs, 2) ); 
wid = full([ sum( coef.* model.SVs ), -model.rho ]);

% 使用ssvm训练的w来测试，看是否线性可分，从而决定是否采用核函数 2015.10.13
if 0
    test_linear_using_ssvm_result( 'fid', positive_sample, negative_sample );
end

%% 提取迁移move样本
positive_sample = [];
negative_sample = [];
num_p_s = 1;
num_n_s = 1;

for t=1:frame-1
    for j=1:n(t)
        if sum( Fij{t}(j,:) )~=0 % 出现迁移则记下正样本
        	k_p = find(Fij{t}(j,:)==1);
            % 样本储存格式
            positive_sample(num_p_s,:) = [ t, j, k_p ];
            num_p_s = num_p_s + 1;
            % 同时选取其他几个中的一个作为负样本
            k_ct = find(Fij{t}(j,:)==0);

            k_n = k_ct( randi(numel(k_ct)) );  % 随机选取一个负样本
            negative_sample(num_n_s,:) = [ t, j, k_n ];
            num_n_s = num_n_s + 1;
        else
            % 没出现迁移则随机选取一个负样本
            k_n = randi(4);       
            negative_sample(num_n_s,:) = [ t, j, k_n ];
            num_n_s = num_n_s + 1;
        end
    end
end
num_p_s = num_p_s -1; % 正样本数目
num_n_s = num_n_s -1; % 负样本数目
    
% 输入到SVM分类器中验证是否线性可分
model = CXSL_SVM( positive_sample, negative_sample, feature_fij, 2/3, 1, 'rand');% 倒数第二位表示测试中负样本数目是正样本的多少倍

coef = repmat( model.sv_coef, 1, size(model.SVs, 2) ); 
wij = full([ sum( coef.* model.SVs ), -model.rho ]);
%
% 速度极慢，以下是运行结果
% ===== Confusion Matrix =====
% 
% 	+1	-1	<-- classified as
% 	651	1	|	+1
% 	6	646	|	-1
% 
% Elapsed time is 107.656502 seconds.
%
% 使用ssvm训练的w来测试，看是否线性可分，从而决定是否采用核函数 2015.10.13
if 0
    test_linear_using_ssvm_result( 'fij', positive_sample, negative_sample );
end

%% 提取appear样本

positive_sample = [];
negative_sample = [];
num_p_s = 1;
num_n_s = 1;

for t=2:frame
    if sum(Fsj{t})~=0 % 如果这一帧中含有正样本
        % 找出正样本
        j_p = find(Fsj{t}(:)==1);
        for i=1:numel(j_p)
            positive_sample(num_p_s,:) = [ t, j_p(i), 1 ];
            num_p_s = num_p_s + 1;
        end
        % 在其他几个中随机抽取相等数目的负样本
        j_ct = setdiff(1:n(t), j_p);
        j_n = j_ct( randi(numel(j_ct), [1 numel(j_p)]) );
        
        for i=1:numel(j_n)
        	negative_sample(num_n_s,:) = [ t, j_n(i), 1 ];
            num_n_s = num_n_s + 1;
        end
    else
        % 随机抽取若干个个负样本    
        j_n = randi(n(t), [1 round(n(t)/20)]);
        for i=1:numel(j_n)
        	negative_sample(num_n_s,:) = [ t, j_n(i), 1 ];
            num_n_s = num_n_s + 1;
        end
    end
end
   
num_p_s = num_p_s -1; % 正样本数目
num_n_s = num_n_s -1; % 负样本数目

model = CXSL_SVM( positive_sample, negative_sample, feature_fsj, 2/3, 2, 'rand');

coef = repmat( model.sv_coef, 1, size(model.SVs, 2) ); 
wsj = full([ sum( coef.* model.SVs ), -model.rho ]);

% 使用ssvm训练的w来测试，看是否线性可分，从而决定是否采用核函数 2015.10.13
if 0
    test_linear_using_ssvm_result( 'fsj', positive_sample, negative_sample );
end

%% 提取disappear样本

positive_sample = [];
negative_sample = [];
num_p_s = 1;
num_n_s = 1;

for t=1:frame-1
    if sum(Fit{t})~=0 % 如果这一帧中含有正样本
        % 找出正样本
        j_p = find(Fit{t}(:)==1);
        for i=1:numel(j_p)
            positive_sample(num_p_s,:) = [ t, j_p(i), 1 ];
            num_p_s = num_p_s + 1;
        end
        % 在其他几个中随机抽取相等数目的负样本
        j_ct = setdiff(1:n(t), j_p);
        j_n = j_ct( randi(numel(j_ct), [1 numel(j_p)]) );
        for i=1:numel(j_n)
        	negative_sample(num_n_s,:) = [ t, j_n(i), 1 ];
            num_n_s = num_n_s + 1;
        end
    else
        % 随机抽取若干个个负样本    
        j_n = randi(n(t), [1 round(n(t)/20)]);
        for i=1:numel(j_n)
        	negative_sample(num_n_s,:) = [ t, j_n(i), 1 ];
            num_n_s = num_n_s + 1;
        end
    end
end

num_p_s = num_p_s -1; % 正样本数目
num_n_s = num_n_s -1; % 负样本数目

model = CXSL_SVM( positive_sample, negative_sample, feature_fit, 2/3, 2, 'rand');

coef = repmat( model.sv_coef, 1, size(model.SVs, 2) ); 
wit = full([ sum( coef.* model.SVs ), -model.rho ]);

% 使用ssvm训练的w来测试，看是否线性可分，从而决定是否采用核函数 2015.10.13
if 0
    test_linear_using_ssvm_result( 'fit', positive_sample, negative_sample );
end

%% 提取merge样本

positive_sample = [];
negative_sample = [];
num_p_s = 1;
num_n_s = 1;

for t=2:frame
    for j=1:n(t)
        if sum( Fmj{t}(j,:) )~=0 % 出现分裂则记下正样本
            mm_p = find(Fmj{t}(j,:)==1);
            % 样本储存格式
            positive_sample(num_p_s,:) = [ t, j, mm_p, candidate_k_last{t}{j,mm_p}];
            num_p_s = num_p_s + 1;
            
            % 同时选取其他几个中的一个作为负样本
            mm_ct = find(Fmj{t}(j,:)==0);  

            mm_n = mm_ct( randi(numel(mm_ct)) );  % 随机选取一个负样本
            negative_sample(num_n_s,:) = [ t, j, mm_n, candidate_k_last{t}{j,mm_n}];
            num_n_s = num_n_s + 1;
            
        else
            % =======================================
            % 设置一个随机过滤器，来减少负样本的采样数量
            if randi(20)~=1
                continue; 
            end % 只有1/10的概率会选择一个负样本
            % =======================================

            mm_n = randi(6);       % 没出现分裂则随机选取一个负样本
            negative_sample(num_n_s,:) = [ t, j, mm_n, candidate_k_last{t}{j,mm_n}];
            num_n_s = num_n_s + 1;
        end
    end
end
num_p_s = num_p_s -1; % 正样本数目
num_n_s = num_n_s -1; % 负样本数目

% 输入到SVM分类器中验证是否线性可分
% 选择 'margin' 模式的时候2/3是失效的，训练样本数目就是边界点数
model = CXSL_SVM( positive_sample, negative_sample, feature_fmj, 2/3, 10, 'rand');% 倒数第二位表示测试中负样本数目是正样本的多少倍

coef = repmat( model.sv_coef, 1, size(model.SVs, 2) ); 
wmj = full([ sum( coef.* model.SVs ), -model.rho ]);
if wmj==0 % 若没有该类型样本，则w会变为0，需要另外设置
    wmj = [zeros(size(feature_fmj{2}{1}))', 0];
end

% 使用ssvm训练的w来测试，看是否线性可分，从而决定是否采用核函数 2015.10.13
if 0
    test_linear_using_ssvm_result( 'fmj', positive_sample, negative_sample );
end

%% 提取split样本

positive_sample = [];
negative_sample = [];
num_p_s = 1;
num_n_s = 1;

for t=1:frame-1
    for j=1:n(t)
        if sum( Fiv{t}(j,:) )~=0 % 出现分裂则记下正样本
            mm_p = find(Fiv{t}(j,:)==1);
            % 样本储存格式
            positive_sample(num_p_s,:) = [ t, j, mm_p, candidate_k_next{t}{j,mm_p}];
            num_p_s = num_p_s + 1;
            
            % 同时选取其他几个中的一个作为负样本
            mm_ct = find(Fiv{t}(j,:)==0);  

            mm_n = mm_ct( randi(numel(mm_ct)) );  % 随机选取一个负样本
            negative_sample(num_n_s,:) = [ t, j, mm_n, candidate_k_next{t}{j,mm_n}];
            num_n_s = num_n_s + 1;
            
        else
            % =======================================
            % 设置一个随机过滤器，来减少负样本的采样数量
            if randi(20)~=1
                continue; 
            end % 只有1/10的概率会选择一个负样本
            % =======================================

            mm_n = randi(6);       % 没出现分裂则随机选取一个负样本
            negative_sample(num_n_s,:) = [ t, j, mm_n, candidate_k_next{t}{j,mm_n}];
            num_n_s = num_n_s + 1;
        end
    end
end
num_p_s = num_p_s -1; % 正样本数目
num_n_s = num_n_s -1; % 负样本数目

% 输入到SVM分类器中验证是否线性可分
% 选择 'margin' 模式的时候2/3是失效的，训练样本数目就是边界点数
model = CXSL_SVM( positive_sample, negative_sample, feature_fiv, 1/2, 10, 'rand');% 倒数第二位表示测试中负样本数目是正样本的多少倍

coef = repmat( model.sv_coef, 1, size(model.SVs, 2) ); 
wiv = full([ sum( coef.* model.SVs ), -model.rho ]);
if wiv==0 % 若没有该类型样本，则w会变为0，需要另外设置
    wiv = [zeros(size(feature_fiv{2}{1}))', 0];
end
% 使用ssvm训练的w来测试，看是否线性可分，从而决定是否采用核函数 2015.10.13
if 0
    test_linear_using_ssvm_result( 'fsj', positive_sample, negative_sample );
end

%% 提取虚景 false detetion 样本（暂时不用）

positive_sample = [];
negative_sample = [];
num_p_s = 1;
num_n_s = 1;

%% 保存各细胞事件的w
if 1
    save([ trackpath, '\结构化学习\initial_w_New.mat'], 'wij', 'wit', 'wid', 'wiv', 'wmj', 'wsj');
end












            
            
        
        
        
