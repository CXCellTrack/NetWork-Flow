%
% ======================================================= %
%
% ʹ��libsvm���������Ƿ����Կɷ�
% �����������ȫ���������в��ԣ��� CXSL_Test_Linear ��������
% ����������Ϻ÷���ĸ��¼� w(����) ������������Ϊ�ṹ��ѧϰ�ĳ�ʼֵ
% 
% ======================================================= %

clear;close all;

dataset = 'training'; % ����ű�һ��ֻ��ѵ������ʹ��
[ segpath trackpath ] = getpath( dataset );

load([ trackpath, '\Pair\Pre_data_New.mat'],'candidate_k_last','candidate_k_next','n','candidate_fij');
% �����׼��GT
load([ trackpath, '\GT\GT_after_hand_tune\GT_Flow_Variables_New.mat']);
% �������������������㣩
load([ trackpath, '\�ṹ��ѧϰ\Feature_New.mat']);
frame = numel(Fmj);    

%% ��ȡ����divide����

positive_sample = [];
negative_sample = [];
num_p_s = 1;
num_n_s = 1;

for t=1:frame-1
    for j=1:n(t)
        if sum( Fid{t}(j,:) )~=0 % ���ַ��������������
            mm_p = find(Fid{t}(j,:)==1);
            % ���������ʽ
            positive_sample(num_p_s,:) = [ t, j, mm_p, candidate_k_next{t}{j,mm_p}];
            num_p_s = num_p_s + 1;
            
            % ͬʱѡȡ���������е�һ����Ϊ������
            mm_ct = find(Fid{t}(j,:)==0);

            mm_n = mm_ct( randi(numel(mm_ct)) );  % ���ѡȡһ��������
            negative_sample(num_n_s,:) = [ t, j, mm_n, candidate_k_next{t}{j,mm_n}];
            num_n_s = num_n_s + 1;
            
        else
            % =======================================
            % ����һ������������������ٸ������Ĳ�������
            if randi(10)~=1
                continue; 
            end % ֻ��1/10�ĸ��ʻ�ѡ��һ��������
            % =======================================

            mm_n = randi(6);       % û���ַ��������ѡȡһ��������
            negative_sample(num_n_s,:) = [ t, j, mm_n, candidate_k_next{t}{j,mm_n}];
            num_n_s = num_n_s + 1;
        end
    end
end
num_p_s = num_p_s -1; % ��������Ŀ
num_n_s = num_n_s -1; % ��������Ŀ

% ʹ��ϵͳ�Դ���svmtrain���з��ࣨ�ٶ������ֵĻ�������ˣ���
[ wid bid ] = CXSL_SVM_new( positive_sample, negative_sample, feature_fid, 1/3);
% bid = bid + 1;
% % ���뵽SVM����������֤�Ƿ����Կɷ�
% % ѡ�� 'margin' ģʽ��ʱ��2/3��ʧЧ�ģ�ѵ��������Ŀ���Ǳ߽����
% [ model bestc bestg ] = CXSL_SVM( positive_sample, negative_sample, feature_fid, 2/3, 3, 'rand');% �����ڶ�λ��ʾ�����и�������Ŀ���������Ķ��ٱ�
% 
% coef = repmat( model.sv_coef, 1, size(model.SVs, 2) ); 
% bid = -model.rho;
% wid = full( sum(coef.* model.SVs) ); % 

% ʹ��ssvmѵ����w�����ԣ����Ƿ����Կɷ֣��Ӷ������Ƿ���ú˺��� 2015.10.13
if 0
    test_linear_using_ssvm_result( 'fid', positive_sample, negative_sample );
end

%% ��ȡǨ��move����
positive_sample = [];
negative_sample = [];
num_p_s = 1;
num_n_s = 1;

for t=1:frame-1
    for j=1:n(t)
        if sum( Fij{t}(j,:) )~=0 % ����Ǩ�������������
        	k_p = find(Fij{t}(j,:)==1);
            % ���������ʽ
            positive_sample(num_p_s,:) = [ t, j, k_p ];
            num_p_s = num_p_s + 1;
            % ͬʱѡȡ���������е�һ����Ϊ������
            k_ct = find(Fij{t}(j,:)==0);

            k_n = k_ct( randi(numel(k_ct)) );  % ���ѡȡһ��������
            negative_sample(num_n_s,:) = [ t, j, k_n ];
            num_n_s = num_n_s + 1;
        else
            % û����Ǩ�������ѡȡһ��������
            k_n = randi(4);       
            negative_sample(num_n_s,:) = [ t, j, k_n ];
            num_n_s = num_n_s + 1;
        end
    end
end
num_p_s = num_p_s -1; % ��������Ŀ
num_n_s = num_n_s -1; % ��������Ŀ
    
[ wij bij ] = CXSL_SVM_new( positive_sample, negative_sample, feature_fij, 1/3);
%
% �ٶȼ��������������н��
% ===== Confusion Matrix =====
% 
% 	+1	-1	<-- classified as
% 	651	1	|	+1
% 	6	646	|	-1
% 
% Elapsed time is 107.656502 seconds.
%
% ʹ��ssvmѵ����w�����ԣ����Ƿ����Կɷ֣��Ӷ������Ƿ���ú˺��� 2015.10.13
if 0
    test_linear_using_ssvm_result( 'fij', positive_sample, negative_sample );
end

%% ��ȡappear����

positive_sample = [];
negative_sample = [];
num_p_s = 1;
num_n_s = 1;

for t=2:frame
    if sum(Fsj{t})~=0 % �����һ֡�к���������
        % �ҳ�������
        j_p = find(Fsj{t}(:)==1);
        for i=1:numel(j_p)
            positive_sample(num_p_s,:) = [ t, j_p(i), 1 ];
            num_p_s = num_p_s + 1;
        end
        % �����������������ȡ�����Ŀ�ĸ�����
        j_ct = setdiff(1:n(t), j_p);
        j_n = j_ct( randi(numel(j_ct), [1 numel(j_p)]) );
        
        for i=1:numel(j_n)
        	negative_sample(num_n_s,:) = [ t, j_n(i), 1 ];
            num_n_s = num_n_s + 1;
        end
    else
        % �����ȡ���ɸ���������    
        j_n = randi(n(t), [1 round(n(t)/20)]);
        for i=1:numel(j_n)
        	negative_sample(num_n_s,:) = [ t, j_n(i), 1 ];
            num_n_s = num_n_s + 1;
        end
    end
end
   
num_p_s = num_p_s -1; % ��������Ŀ
num_n_s = num_n_s -1; % ��������Ŀ

[ wsj bsj ] = CXSL_SVM_new( positive_sample, negative_sample, feature_fsj, 1/3);
% bsj = bsj + 1;
% ʹ��ssvmѵ����w�����ԣ����Ƿ����Կɷ֣��Ӷ������Ƿ���ú˺��� 2015.10.13
if 0
    test_linear_using_ssvm_result( 'fsj', positive_sample, negative_sample );
end

%% ��ȡdisappear����

positive_sample = [];
negative_sample = [];
num_p_s = 1;
num_n_s = 1;

for t=1:frame-1
    if sum(Fit{t})~=0 % �����һ֡�к���������
        % �ҳ�������
        j_p = find(Fit{t}(:)==1);
        for i=1:numel(j_p)
            positive_sample(num_p_s,:) = [ t, j_p(i), 1 ];
            num_p_s = num_p_s + 1;
        end
        % �����������������ȡ�����Ŀ�ĸ�����
        j_ct = setdiff(1:n(t), j_p);
        j_n = j_ct( randi(numel(j_ct), [1 numel(j_p)]) );
        for i=1:numel(j_n)
        	negative_sample(num_n_s,:) = [ t, j_n(i), 1 ];
            num_n_s = num_n_s + 1;
        end
    else
        % �����ȡ���ɸ���������    
        j_n = randi(n(t), [1 round(n(t)/20)]);
        for i=1:numel(j_n)
        	negative_sample(num_n_s,:) = [ t, j_n(i), 1 ];
            num_n_s = num_n_s + 1;
        end
    end
end

num_p_s = num_p_s -1; % ��������Ŀ
num_n_s = num_n_s -1; % ��������Ŀ

[ wit bit ] = CXSL_SVM_new( positive_sample, negative_sample, feature_fit, 1/3);
% bit = bit + 1;
% ʹ��ssvmѵ����w�����ԣ����Ƿ����Կɷ֣��Ӷ������Ƿ���ú˺��� 2015.10.13
if 0
    test_linear_using_ssvm_result( 'fit', positive_sample, negative_sample );
end

%% ��ȡmerge����

positive_sample = [];
negative_sample = [];
num_p_s = 1;
num_n_s = 1;

for t=2:frame
    for j=1:n(t)
        if sum( Fmj{t}(j,:) )~=0 % ���ַ��������������
            mm_p = find(Fmj{t}(j,:)==1);
            % ���������ʽ
            positive_sample(num_p_s,:) = [ t, j, mm_p, candidate_k_last{t}{j,mm_p}];
            num_p_s = num_p_s + 1;
            
            % ͬʱѡȡ���������е�һ����Ϊ������
            mm_ct = find(Fmj{t}(j,:)==0);  

            mm_n = mm_ct( randi(numel(mm_ct)) );  % ���ѡȡһ��������
            negative_sample(num_n_s,:) = [ t, j, mm_n, candidate_k_last{t}{j,mm_n}];
            num_n_s = num_n_s + 1;
            
        else
            % =======================================
            % ����һ������������������ٸ������Ĳ�������
            if randi(20)~=1
                continue; 
            end % ֻ��1/10�ĸ��ʻ�ѡ��һ��������
            % =======================================

            mm_n = randi(6);       % û���ַ��������ѡȡһ��������
            negative_sample(num_n_s,:) = [ t, j, mm_n, candidate_k_last{t}{j,mm_n}];
            num_n_s = num_n_s + 1;
        end
    end
end
num_p_s = num_p_s -1; % ��������Ŀ
num_n_s = num_n_s -1; % ��������Ŀ

[ wmj bmj ] = CXSL_SVM_new( positive_sample, negative_sample, feature_fmj, 1/3);
% bmj = bmj + 1;

% ʹ��ssvmѵ����w�����ԣ����Ƿ����Կɷ֣��Ӷ������Ƿ���ú˺��� 2015.10.13
if 0
    test_linear_using_ssvm_result( 'fmj', positive_sample, negative_sample );
end

%% ��ȡsplit����

positive_sample = [];
negative_sample = [];
num_p_s = 1;
num_n_s = 1;

for t=1:frame-1
    for j=1:n(t)
        if sum( Fiv{t}(j,:) )~=0 % ���ַ��������������
            mm_p = find(Fiv{t}(j,:)==1);
            % ���������ʽ
            positive_sample(num_p_s,:) = [ t, j, mm_p, candidate_k_next{t}{j,mm_p}];
            num_p_s = num_p_s + 1;
            
            % ͬʱѡȡ���������е�һ����Ϊ������
            mm_ct = find(Fiv{t}(j,:)==0);  

            mm_n = mm_ct( randi(numel(mm_ct)) );  % ���ѡȡһ��������
            negative_sample(num_n_s,:) = [ t, j, mm_n, candidate_k_next{t}{j,mm_n}];
            num_n_s = num_n_s + 1;
            
        else
            % =======================================
            % ����һ������������������ٸ������Ĳ�������
            if randi(20)~=1
                continue; 
            end % ֻ��1/10�ĸ��ʻ�ѡ��һ��������
            % =======================================

            mm_n = randi(6);       % û���ַ��������ѡȡһ��������
            negative_sample(num_n_s,:) = [ t, j, mm_n, candidate_k_next{t}{j,mm_n}];
            num_n_s = num_n_s + 1;
        end
    end
end
num_p_s = num_p_s -1; % ��������Ŀ
num_n_s = num_n_s -1; % ��������Ŀ

[ wiv biv ] = CXSL_SVM_new( positive_sample, negative_sample, feature_fiv, 1/3);

% ʹ��ssvmѵ����w�����ԣ����Ƿ����Կɷ֣��Ӷ������Ƿ���ú˺��� 2015.10.13
if 0
    test_linear_using_ssvm_result( 'fsj', positive_sample, negative_sample );
end

%% �����ϸ���¼���w
if 1
    save([ trackpath, '\�ṹ��ѧϰ\initial_w_New.mat'],...
        'wij', 'wit', 'wid', 'wiv', 'wmj', 'wsj',...
        'bij', 'bit', 'bid', 'biv', 'bmj', 'bsj');
end












            
            
        
        
        
