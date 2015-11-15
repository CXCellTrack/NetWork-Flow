function  [ model bestc bestg ] = CXSL_SVM( positive_sample, negative_sample, feature, n_cut, num_n_p, option)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%
%   n_cut   ����������ѵ�������ı���
%   num_n_p ָ��������Ŀ���������Ķ��ٱ����ڲ�����ʹ���������
%   option  Ŀǰ��2��ѡ�� 'mindist'   Ϊ��ȡ���������о����������Ϊѵ����
%                          'rand'      Ϊ�������
%                          'margin'    Ϊ��ȡ���������ֲ�������ı߽����Ϊѵ����
tic;

num_p_s = size(positive_sample, 1); % ��������Ŀ
num_n_s = size(negative_sample, 1); % ��������Ŀ

% ����ѵ���������� n_cut;
num_train_p = round(num_p_s * n_cut);
num_train_n = num_train_p;
% ѵ����������ѡ����ͬ��Ŀ������ѡ������Ϊ���������ɱ�
num_test_p = num_p_s - num_train_p;
num_test_n = num_test_p * num_n_p;

if strcmp(option, 'mindist')
%% 1.������������Ѱ�������ƽ�����Ϊѵ������Ч���

p_feature = cell(1,1);
for i=1:num_p_s
    e_info = positive_sample(i,:);
    p_feature{i,1} = feature{e_info(1)}{e_info(2), e_info(3)}';
end
p_feature = cell2mat(p_feature); % ���ȫ������������

n_feature = cell(1,1);
for i=1:num_n_s
    e_info = negative_sample(i,:);
    n_feature{i,1} = feature{e_info(1)}{e_info(2), e_info(3)}';
end
n_feature = cell2mat(n_feature)'; % ���ȫ�帺��������

distance = dist(p_feature, n_feature); % ����������֮���ŷ�Ͼ���

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

[~, ind_p] = sortrows(p_with_n, 1);% �������С��������
[~, ind_n] = sortrows(n_with_p, 1);% �������С��������

ind_p_train = ind_p(1:num_train_p);
ind_n_train = ind_n(1:num_train_n);
train_label = [ ones(num_train_p,1);-ones(num_train_n,1) ];
train_feature = [ p_feature(ind_p_train,:); n_feature(:,ind_n_train)' ];

ind_p_test = setdiff(1:num_p_s, ind_p_train)';
ind_n_test = setdiff(1:num_n_s, ind_n_train)';
% �����Ը���������˳�򣬴��г�ȡ num_test_n ��
ind_n_test = ind_n_test( randperm(numel(ind_n_test)) );
ind_n_test = ind_n_test(1:num_test_n);
test_label = [ ones(num_test_p,1);-ones(num_test_n,1) ];

test_feature = [ p_feature(ind_p_test,:); n_feature(:,ind_n_test)' ];

elseif strcmp(option, 'margin')
%% 2.�������Ե������(����Ҫ����ѵ�����Ա��� n_cut )(Ч������)

p_feature = cell(1,1);
for i=1:num_p_s
    e_info = positive_sample(i,:);
    p_feature{i,1} = feature{e_info(1)}{e_info(2), e_info(3)}';
end
p_feature = cell2mat(p_feature); % ���ȫ������������

n_feature = cell(1,1);
for i=1:num_n_s
    e_info = negative_sample(i,:);
    n_feature{i,1} = feature{e_info(1)}{e_info(2), e_info(3)}';
end
n_feature = cell2mat(n_feature); % ���ȫ�帺��������

% �������Ļ�����
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
ind_p_train = unique(ind_maxmin_p); % �ҳ��������ı߽�
ind_n_train = unique(ind_maxmin_n); % �ҳ��������ı߽�

num_train_p = numel(ind_p_train); % ȡ�߽��ϵĵ����ѵ��
num_train_n = numel(ind_n_train);
% ����ѵ�������ͱ�ǩ
train_label = [ ones(num_train_p,1);-ones(num_train_n,1) ];
train_feature = [ p_feature(ind_p_train,:); n_feature(ind_n_train,:) ];

ind_p_test = setdiff(1:num_p_s, ind_p_train)';
ind_n_test = setdiff(1:num_n_s, ind_n_train)';
% �����Ը���������˳�򣬴��г�ȡ num_test_n ��
num_test_p = numel(ind_p_test);
num_test_n = num_test_p * num_n_p;

ind_n_test = ind_n_test( randperm(numel(ind_n_test)) );
ind_n_test = ind_n_test(1:num_test_n);

test_label = [ ones(num_test_p,1);-ones(num_test_n,1) ];
test_feature = [ p_feature(ind_p_test,:); n_feature(ind_n_test,:) ];

elseif strcmp(option, 'rand')
%% 3.��ͨ���������
% ��������������˳��
rand_indp = randperm(num_p_s);
rand_indn = randperm(num_n_s);
% ��ȡ��Ŀ��ȵ�������������ѵ��
ind_train_p = rand_indp(1:num_train_p);
ind_train_n = rand_indn(1:num_train_n);
% ��������ʣ�µĵ�ȫ����Ϊ����������
ind_test_p = mysetdiff(1:num_p_s, ind_train_p);
% ��������ʣ�µĵ�����һ�������Ĳ��Ը�����
all_test_n = mysetdiff(1:num_n_s, ind_train_n);
ind_test_n = all_test_n(1:num_test_n);

% ��������ѵ����������������
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
% ��������������������������
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
    error('�������������� rand/margin/mindist �е�һ��');
end
%% ��ʼѵ��SVM������
disp('start svm training...');
if 0
    % Ѱ�����ŵ�cg
    [~, bestc, bestg] = SVMcgForClass(train_label,train_feature,-10,10,-10,10,3,1,1);
else
    bestc = 1;
    bestg = 1;
end

model = svmtrain(train_label, train_feature, ['-t 0 -c ',num2str(bestc),' -g ',num2str(bestg)]); % �� -b 1 �ɵõ��������
disp('start svm predicting...');
[ predicted_label ] = svmpredict(test_label, test_feature, model); % �� -b 1 �ɵõ��������

%% �����������
cm = zeros(2,2);
tp = find(test_label(1:num_test_p) == predicted_label(1:num_test_p));
cm(1,1) = numel(tp);
cm(2,1) = num_test_p - cm(1,1);
tp = find(test_label(num_test_p+1:end) == predicted_label(num_test_p+1:end));
cm(2,2) = numel(tp);
cm(1,2) = num_test_n - cm(2,2);

% ��ӡ��Ϣ
fprintf('\n===== Confusion Matrix =====\n\n');
fprintf('\t\t\t\t��ʵ���\n');
fprintf('\t\t\t\t+1\t-1\n');
fprintf('\t������+1\t%d\t%d\n',cm(1,:));
fprintf('\t������-1\t%d\t%d\n\n',cm(2,:));

precision = cm(1,1)/sum(cm(1,:));
recall = cm(1,1)/sum(cm(:,1));
F_measure = 2*precision*recall/(precision+recall);
fprintf('\tprecision = %0.3f\n', precision);
fprintf('\trecall    = %0.3f\n', recall);
fprintf('\tF_measure = %0.3f\n', F_measure);
toc;

end

