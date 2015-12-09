function test_linear_using_ssvm_result( event, positive_sample, negative_sample )

% 这个函数：使用ssvm得到的w来分割各事件，如果能分开，则说明特征线性可分；
% 如果不能分开，则说明特征线性可分性差，如果使用核可能有较大的提升
% 2015.10.13

% 1、载入SSVM训练得到的最佳 w        
[ ~, traintrackpath ] = getpath( 'training' );
if 0
    disp('  载入之前 SSVM 训练得到的最佳 w...');
    load([ traintrackpath, '\结构化学习\SSVM_Best_W_New.mat']);
else % 想要复现 excel 中记载的以前的实验结果，只需要手动填写 w_best 即可
    disp('  载入手动填写的w...');
    w_best = [-7.935998125	1.720791701	0.18538702	0.054008597	-0.108097034	-0.968294343	0.636409321	0.646062796	-7.370071068	-2.567381363	-0.362143169	0.18350773	-0.813164628	1.280271937	-2.600506538	7.796467214	0.54077495	-0.118970058	-2.03110677	-2.38126945	0.615545764	0.114766087	1.680866361	-10.65127626	4.482456247	2.009258903	-0.642628466	-1.814612604	-3.308992643	1.293138786	0.87139507	1.664626957	-0.014863296	-0.845479662	-3.942422168	-0.481027921	-1.496917921	0.163797962	-0.31106576	-5.608061028]';
end

%% 2、从w_best中将各事件的w分割出来
load([ traintrackpath, '\结构化学习\initial_w_New.mat']);
df = zeros(6,1);
df(1) = 1;
df(2) = numel(wij)+1;
df(3) = numel(wij)+numel(wit)+1;
df(4) = numel(wij)+numel(wit)+numel(wid)+1;
df(5) = numel(wij)+numel(wit)+numel(wid)+numel(wiv)+1;
df(6) = numel(wij)+numel(wit)+numel(wid)+numel(wiv)+numel(wmj)+1;

wij = w_best( df(1):df(2)-1 )'; % 动态分割赋值，以后修改特征就无需修改这里了
wit = w_best( df(2):df(3)-1 )';
wid = w_best( df(3):df(4)-1 )';
wiv = w_best( df(4):df(5)-1 )';
wmj = w_best( df(5):df(6)-1 )';
wsj = w_best( df(6):end )';

%% 3、根据输入字符选择不同的事件
load([ traintrackpath, '\结构化学习\Feature_Plus_New.mat']);
switch event
    case 'fij'
        feature = feature_fij_p;
        w = wij;
    case 'fit'
        feature = feature_fit_p;
        w = wit;
    case 'fid'
        feature = feature_fid_p;
        w = wid;
    case 'fiv'
        feature = feature_fiv_p;
        w = wiv;
    case 'fmj'
        feature = feature_fmj_p;
        w = wmj;
    case 'fsj'
        feature = feature_fsj_p;
        w = wsj;
end

%% 制造正负训练样本的特征矩阵
num_p = size(positive_sample,1);
num_n = size(negative_sample,1);

label = [ ones(num_p,1);-ones(num_n,1) ];
sample_feature = cell(num_p + num_n,1);
for i=1:num_p
    e_info = positive_sample(i,:);
    sample_feature{i,1} = feature{e_info(1)}{e_info(2), e_info(3)}';
end
for i=1:num_n
    e_info = negative_sample(i,:);
    sample_feature{i+num_p,1} = feature{e_info(1)}{e_info(2), e_info(3)}';
end
sample_feature = cell2mat(sample_feature);

%% 5、使用<w,feature>得到分类结果
predicted_label = (w*sample_feature')';
predicted_label(predicted_label>=0) = 1;
predicted_label(predicted_label<0) = -1;

cm = zeros(2,2);
tp = find(label(1:num_p)==predicted_label(1:num_p));
cm(1,1) = numel(tp);
cm(2,1) = num_p - cm(1,1);
tp = find(label(num_p+1:end)==predicted_label(num_p+1:end));
cm(2,2) = numel(tp);
cm(1,2) = num_n - cm(2,2);

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












