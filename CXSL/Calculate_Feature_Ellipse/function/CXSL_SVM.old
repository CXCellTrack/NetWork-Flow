function  [ model bestc bestg ] = CXSL_SVM( positive_sample, negative_sample, feature, n_cut, num_n_p, option)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%
%   n_cut   是正样本中训练样本的比例
%   num_n_p 指负样本数目是正样本的多少倍，在测试中使用这个倍数
%   option  目前有2种选择： 'mindist'   为抽取正负样本中距离最近的作为训练集
%                          'rand'      为随机抽样
%                          'margin'    为抽取正负样本分布超球体的边界点作为训练集
tic;

num_p_s = size(positive_sample, 1); % 正样本数目
num_n_s = size(negative_sample, 1); % 负样本数目

% 设置训练样本比例 n_cut;
num_train_p = round(num_p_s * n_cut);
num_train_n = num_train_p;
% 训练正负样本选择相同数目，测试选择负样本为正样本若干倍
num_test_p = num_p_s - num_train_p;
num_test_n = num_test_p * num_n_p;

if strcmp(option, 'mindist')
%% 1.按照特征距离寻找最靠近超平面的作为训练集（效果差）

p_feature = cell(1,1);
for i=1:num_p_s
    e_info = positive_sample(i,:);
    p_feature{i,1} = feature{e_info(1)}{e_info(2), e_info(3)}';
end
p_feature = cell2mat(p_feature); % 求出全体正样本特征

n_feature = cell(1,1);
for i=1:num_n_s
    e_info = negative_sample(i,:);
    n_feature{i,1} = feature{e_info(1)}{e_info(2), e_info(3)}';
end
n_feature = cell2mat(n_feature)'; % 求出全体负样本特征

distance = dist(p_feature, n_feature); % 求正负样本之间的欧氏距离

p_with_n = zeros(num_p_s,1);
for h=1:num_p_s
    p_with_n(h,1) = min(distance(h,:));
    p_with_n(h,2) = find( distance(h,:) == p_with_n(h,1) );
end
n_with_p = zeros(num_n_s,1);
for w=1:num_n_s
    n_with_p(w,1) = min(distance(:,w));
    n_with_p(w,2) = find( distance(:,w) == n_with_p(w,1) );
end

[~, ind_p] = sortrows(p_with_n, 1);% 按距离从小到大排列
[~, ind_n] = sortrows(n_with_p, 1);% 按距离从小到大排列

ind_p_train = ind_p(1:num_train_p);
ind_n_train = ind_n(1:num_train_n);
train_label = [ ones(num_train_p,1);-ones(num_train_n,1) ];
train_feature = [ p_feature(ind_p_train,:); n_feature(:,ind_n_train)' ];

ind_p_test = setdiff(1:num_p_s, ind_p_train)';
ind_n_test = setdiff(1:num_n_s, ind_n_train)';
% 将测试负样本打乱顺序，从中抽取 num_test_n 个
ind_n_test = ind_n_test( randperm(numel(ind_n_test)) );
ind_n_test = ind_n_test(1:num_test_n);
test_label = [ ones(num_test_p,1);-ones(num_test_n,1) ];

test_feature = [ p_feature(ind_p_test,:); n_feature(:,ind_n_test)' ];

elseif strcmp(option, 'margin')
%% 2.超球体边缘抽样法(不需要分配训练测试比例 n_cut )(效果还行)

p_feature = cell(1,1);
for i=1:num_p_s
    e_info = positive_sample(i,:);
    p_feature{i,1} = feature{e_info(1)}{e_info(2), e_info(3)}';
end
p_feature = cell2mat(p_feature); % 求出全体正样本特征

n_feature = cell(1,1);
for i=1:num_n_s
    e_info = negative_sample(i,:);
    n_feature{i,1} = feature{e_info(1)}{e_info(2), e_info(3)}';
end
n_feature = cell2mat(n_feature); % 求出全体负样本特征

% 特征中心化处理
p_feature_o = zscore(p_feature);
n_feature_o = zscore(n_feature);
ind_maxmin_p = zeros(16,1);
ind_maxmin_n = zeros(16,1);
for i=1:8
    tp = p_feature_o(:,i);
    ind_maxmin_p(i) = find(tp == max(tp), 1);
    ind_maxmin_p(i+8) = find(tp == min(tp), 1);
    
    tp = n_feature_o(:,i);
    ind_maxmin_n(i) = find(tp == max(tp), 1);
    ind_maxmin_n(i+8) = find(tp == min(tp), 1);
end
ind_p_train = unique(ind_maxmin_p); % 找出正样本的边界
ind_n_train = unique(ind_maxmin_n); % 找出负样本的边界

num_train_p = numel(ind_p_train); % 取边界上的点进行训练
num_train_n = numel(ind_n_train);
% 制造训练样本和标签
train_label = [ ones(num_train_p,1);-ones(num_train_n,1) ];
train_feature = [ p_feature(ind_p_train,:); n_feature(ind_n_train,:) ];

ind_p_test = setdiff(1:num_p_s, ind_p_train)';
ind_n_test = setdiff(1:num_n_s, ind_n_train)';
% 将测试负样本打乱顺序，从中抽取 num_test_n 个
num_test_p = numel(ind_p_test);
num_test_n = num_test_p * num_n_p;

ind_n_test = ind_n_test( randperm(numel(ind_n_test)) );
ind_n_test = ind_n_test(1:num_test_n);

test_label = [ ones(num_test_p,1);-ones(num_test_n,1) ];
test_feature = [ p_feature(ind_p_test,:); n_feature(ind_n_test,:) ];

elseif strcmp(option, 'rand')
%% 3.普通随机抽样法
% 打乱正负样本的顺序
rand_indp = randperm(num_p_s);
rand_indn = randperm(num_n_s);
% 抽取数目相等的正负样本用作训练
ind_train_p = rand_indp(1:num_train_p);
ind_train_n = rand_indn(1:num_train_n);
% 正样本中剩下的的全部作为测试正样本
ind_test_p = mysetdiff(1:num_p_s, ind_train_p);
% 负样本中剩下的的生成一定数量的测试负样本
all_test_n = mysetdiff(1:num_n_s, ind_train_n);
ind_test_n = all_test_n(1:num_test_n);

% 制造正负训练样本的特征矩阵
train_label = [ ones(num_train_p,1);-ones(num_train_n,1) ];

train_feature = cell(num_train_p + num_train_n,1);
for i=1:num_train_p
    e_info = positive_sample(ind_train_p(i),:);
    train_feature{i,1} = feature{e_info(1)}{e_info(2), e_info(3)}';
end
for i=1:num_train_n
    e_info = negative_sample(ind_train_n(i),:);
    train_feature{i+num_train_p,1} = feature{e_info(1)}{e_info(2), e_info(3)}';
end
train_feature = cell2mat(train_feature);
% 制造正负测试样本的特征矩阵
test_label = [ ones(num_test_p,1);-ones(num_test_n,1) ];

test_feature = cell(num_test_p + num_test_n,1);
for i=1:num_test_p
    e_info = positive_sample(ind_test_p(i),:);
    test_feature{i,1} = feature{e_info(1)}{e_info(2), e_info(3)}';
end
for i=1:num_test_n
    e_info = negative_sample(ind_test_n(i),:);
    test_feature{i+num_test_p,1} = feature{e_info(1)}{e_info(2), e_info(3)}';
end
test_feature = cell2mat(test_feature);

else
    error('采样方法必须是 rand/margin/mindist 中的一种');
end
%% 开始训练SVM并测试
disp('start svm training...');
if 0
    % 寻找最优的cg
    [~, bestc, bestg] = SVMcgForClass(train_label,train_feature,-10,10,-10,10,3,1,1);
else
    bestc = 1;
    bestg = 1;
end

model = svmtrain(train_label, train_feature, ['-t 0 -c ',num2str(bestc),' -g ',num2str(bestg)]); % 加 -b 1 可得到概率输出
disp('start svm predicting...');
[ predicted_label ] = svmpredict(test_label, test_feature, model); % 加 -b 1 可得到概率输出

%% 输出混淆矩阵
cm = zeros(2,2);
tp = find(test_label(1:num_test_p) == predicted_label(1:num_test_p));
cm(1,1) = numel(tp);
cm(2,1) = num_test_p - cm(1,1);
tp = find(test_label(num_test_p+1:end) == predicted_label(num_test_p+1:end));
cm(2,2) = numel(tp);
cm(1,2) = num_test_n - cm(2,2);

% 打印信息
fprintf('\n===== Confusion Matrix =====\n\n');
fprintf('\t\t\t\t真实结果\n');
fprintf('\t\t\t\t+1\t-1\n');
fprintf('\t分类结果+1\t%d\t%d\n',cm(1,:));
fprintf('\t分类结果-1\t%d\t%d\n\n',cm(2,:));

precision = cm(1,1)/sum(cm(1,:));
recall = cm(1,1)/sum(cm(:,1));
F_measure = 2*precision*recall/(precision+recall);
fprintf('\tprecision = %0.3f\n', precision);
fprintf('\trecall    = %0.3f\n', recall);
fprintf('\tF_measure = %0.3f\n', F_measure);
toc;

end

