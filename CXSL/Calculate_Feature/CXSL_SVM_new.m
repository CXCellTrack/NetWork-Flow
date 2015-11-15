function  [ W B ] = CXSL_SVM_new( positive_sample, negative_sample, feature, n_cut)

num_p_s = size(positive_sample, 1); % 正样本数目
num_n_s = size(negative_sample, 1); % 负样本数目

all_feature = cell(num_p_s + num_n_s,1);
for i=1:num_p_s
    e_info = positive_sample(i,:);
    all_feature{i,1} = feature{e_info(1)}{e_info(2), e_info(3)}';
end
for i=1:num_n_s
    e_info = negative_sample(i,:);
    all_feature{i+num_n_s,1} = feature{e_info(1)}{e_info(2), e_info(3)}';
end
all_feature = cell2mat(all_feature);
all_label = [ones(num_p_s,1); -ones(num_n_s,1)];

P = cvpartition(all_label, 'holdout', n_cut);

if 0
    % 使用系统自带的svmtrain
%     svmStruct = svmtrain(all_feature(P.training,:), all_label(P.training), 'autoscale', false);
%     Predicted_label = svmclassify(svmStruct, all_feature(P.test,:));
%     conMat = confusionmat(all_label(P.test), Predicted_label);
% 
%     % 输出w和b，注意这里是w*x+b，与libsvm不太一样
%     W = svmStruct.Alpha'*svmStruct.SupportVectors;
%     B = svmStruct.Bias;
else
    % 使用libsvm
    disp('start svm training...');
    if 0
        % 寻找最优的cg
        [~, bestc, bestg] = SVMcgForClass(train_label,train_feature,-10,10,-10,10,3,1,1);
    else
        bestc = 1;
        bestg = 1;
    end
   
    model = svmtrain(all_label(P.training), all_feature(P.training,:), ['-t 0 -c ',num2str(bestc),' -g ',num2str(bestg)]); % 加 -b 1 可得到概率输出
    disp('start svm predicting...');
%     model.rho = model.rho - 1;
    [ predicted_label ] = svmpredict(all_label(P.test), all_feature(P.test,:), model); % 加 -b 1 可得到概率输出
    conMat = confusionmat(all_label(P.test), predicted_label);
    
    W = full(model.sv_coef'*model.SVs);
    B = -model.rho;
end


%% 将混淆矩阵变为标准形式
TP = conMat(2,2);
conMat(2,2) = conMat(1,1);
conMat(1,1) = TP;
% 打印信息
fprintf('\n===== Confusion Matrix =====\n\n');
fprintf('\t\t\t\t真实结果\n');
fprintf('\t\t\t\t+1\t-1\n');
fprintf('\t分类结果+1\t%d\t%d\n',conMat(1,:));
fprintf('\t分类结果-1\t%d\t%d\n\n',conMat(2,:));

precision = conMat(1,1)/sum(conMat(1,:));
recall = conMat(1,1)/sum(conMat(:,1));
F_measure = 2*precision*recall/(precision+recall);
fprintf('\tprecision = %0.3f\n', precision);
fprintf('\trecall    = %0.3f\n', recall);
fprintf('\tF_measure = %0.3f\n', F_measure);






