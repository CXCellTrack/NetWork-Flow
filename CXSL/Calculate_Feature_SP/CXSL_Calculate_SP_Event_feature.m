% =========================================================================
% 2015.12.9
% 
% 计算 Superpixel 假说各个事件的特征
%
% =========================================================================

clear;close all
if 0
    dataset = 'competition';
else
    dataset = 'training';
end
[ segpath trackpath ] = getpath( dataset );

pre_data_addr = [ trackpath, '\Pair\Pre_data_New.mat'];
load( pre_data_addr );

% 原始维度特征
mkdir([trackpath, '\结构化学习']);
Feature_New_addr = [ trackpath, '\结构化学习\Feature_New.mat'];
% 增广特征
Feature_Plus_New_addr = [ trackpath, '\结构化学习\Feature_Plus_New.mat']; 

%% 计算单个SP的特征，加入到自身的struct中
disp('  计算超像素自身特征...');
tic;
SPF = CXSL_Calculate_SP_feature( SuperPixel, segpath, n ); % SPF为Superpixel的feature
toc;

frame = numel(SPF);
frame = 10;

%% 下面开始计算前后2帧椭圆特征差异
disp('  计算各细胞事件特征...');

%% ================== move ================== 迁移特征全是特征差的形式
%
% 一共12维特征
%
%###########################################  时间花费 32秒 
tic;
feature_fij = cell(frame-1,1);
for t=1:frame-1
    for j=1:n(t)
        spj = SPF{t}{j};
        
        for mm=1:4
            % 映射寻找 k
            k = candidate_fij{t}(j,mm);
            % 候选假说spk
            spk = SPF{t+1}{k};
            % =========================================================== %
            diff_size = abs(spk.Area - spj.Area); % 面积差的绝对值
            distance = norm(spk.Centroid - spj.Centroid); % 位置差异使用欧氏距离
            % 一系列标量特征使用差绝对值
            somej = struct2cell(spj);
            somek = struct2cell(spk);
            diff_some = abs(cell2mat(somej(3:11)) - cell2mat(somek(3:11)));
            % 计算灰度直方图的EMD
            hist_j = spj.HistIntensity;
            hist_k = spk.HistIntensity;
            % --- 调用外部函数 emd_hat_gd_metric_mex 计算 EM 距离
            EMD = emd_hat_gd_metric_mex(hist_j, hist_k, ones(numel(hist_j), numel(hist_j)) ,0 ,1);

            % 合成move特征，最后补上一个1，作为增广特征 2015.6.24（后来移到mapmimax归一化中进行增广的处理）
            feature_fij{t}{j,mm} = [ diff_size; distance; diff_some; EMD ];

        end      
    end
end
toc;

%% ==================  divide ================== 
%
% 一共24维特征
%
%########################################### 
tic;
feature_fid = cell(frame-1,1);
for t=1:frame-1
    for j=1:n(t)
        father = SPF{t}{j};
        
        for mm=1:6
            sons = candidate_k_next{t}{j,mm};
            son1 = SPF{t+1}{sons(1)};
            son2 = SPF{t+1}{sons(2)};
            
            % 灰度和差异
            diff_intensity_sum = father.SumIntensity - son1.SumIntensity - son2.SumIntensity;
            % 面积和差异
            diff_size_sum = father.Area - son1.Area - son2.Area;
            % 距离和
            distance = norm(father.Centroid - son1.Centroid) + norm(father.Centroid - son2.Centroid);
            % 角度，使用向量内积来求角度
            v1 = son1.Centroid - father.Centroid;
            v2 = son2.Centroid - father.Centroid;
            angle_patten = acosd(dot(v1,v2)/( norm(v1)*norm(v2) ));
            
            % 2个子细胞的差异（一系列标量）
            some_son1 = struct2cell(son1);
            some_son2 = struct2cell(son2);
            diff_son = abs(cell2mat(some_son1(3:11)) - cell2mat(some_son2(3:11)));
            
            % 父的部分信息
            some_father = struct2cell(father); % 除了2个向量特征其他都加进来
            father_fe = cell2mat(some_father([1,3:11,13]));
            
            % 合成分裂特征，最后补上一个1，作为增广特征 2015.6.24
            feature_fid{t}{j,mm} = [ diff_intensity_sum; diff_size_sum; distance; angle_patten;...
                diff_son; father_fe ];
        end
    end
end
toc;   

%% ==================  appear ================== 
%
% 一共11维特征
%
%########################################### 
tic
feature_fsj = cell(frame-1,1);
for t=2:frame
    for j=1:n(t)
        spj = SPF{t}{j};
        some = struct2cell(spj); % 除了2个向量特征其他都加进来
        feature = cell2mat(some([1,3:11,13]));
        
        % 最后补上一个1，作为增广特征 2015.6.24
        feature_fsj{t}{j,1} = feature;
    end
end
toc

%% ==================  disappear ================== 
%
% 一共11维特征
%
%########################################### 
tic
feature_fit = cell(frame-1,1);
for t=1:frame-1
    for j=1:n(t)
        spj = SPF{t}{j};
        some = struct2cell(spj); % 除了2个向量特征其他都加进来
        feature = cell2mat(some([1,3:11,13]));
        
         % 最后补上一个1，作为增广特征 2015.6.24
        feature_fit{t}{j,1} = feature;
    end
end
toc

%% ==================  split ================== 
%
% 一共23维特征
%
%########################################### 
tic
feature_fiv = cell(frame-1,1);
for t=1:frame-1
    for j=1:n(t)
        father = SPF{t}{j};
        
        for mm=1:6
            sons = candidate_k_next{t}{j,mm};
            son1 = SPF{t+1}{sons(1)};
            son2 = SPF{t+1}{sons(2)};
            
            % 灰度和差异（这里有abs，与divide不同）
            diff_intensity_sum = abs(father.SumIntensity - son1.SumIntensity - son2.SumIntensity);
            % 面积和差异（这里有abs，与divide不同）
            diff_size_sum = abs(father.Area - son1.Area - son2.Area);
            % 距离和
            distance = norm(father.Centroid - son1.Centroid) + norm(father.Centroid - son2.Centroid);
            
            % 2个子细胞的差异（一系列标量）
            some_son1 = struct2cell(son1);
            some_son2 = struct2cell(son2);
            diff_son = abs(cell2mat(some_son1(3:11)) - cell2mat(some_son2(3:11)));
            
            % 父的部分信息
            some_father = struct2cell(father); % 除了2个向量特征其他都加进来
            father_fe = cell2mat(some_father([1,3:11,13]));
            
            % 合成split特征，最后补上一个1，作为增广特征 2015.6.24
            feature_fiv{t}{j,mm} = [ diff_intensity_sum; diff_size_sum; distance;...
                diff_son; father_fe ];
        end
    end
end
toc

%% ==================  merge ================== 
%
% 一共23维特征
%
%########################################### 
tic
feature_fmj = cell(frame-1,1);
for t=2:frame
    for j=1:n(t)
        father = SPF{t}{j};
        
        for mm=1:6
            sources = candidate_k_last{t}{j,mm};
            son1 = SPF{t-1}{sources(1)};
            son2 = SPF{t-1}{sources(2)};
            
            % 灰度和差异（这里有abs，与divide不同）
            diff_intensity_sum = abs(father.SumIntensity - son1.SumIntensity - son2.SumIntensity);
            % 面积和差异（这里有abs，与divide不同）
            diff_size_sum = abs(father.Area - son1.Area - son2.Area);
            % 距离和
            distance = norm(father.Centroid - son1.Centroid) + norm(father.Centroid - son2.Centroid);
            
            % 2个子细胞的差异（一系列标量）
            some_son1 = struct2cell(son1);
            some_son2 = struct2cell(son2);
            diff_son = abs(cell2mat(some_son1(3:11)) - cell2mat(some_son2(3:11)));
            
            % 父的部分信息
            some_father = struct2cell(father); % 除了2个向量特征其他都加进来
            father_fe = cell2mat(some_father([1,3:11,13]));
            
            % 合成merge特征，最后补上一个1，作为增广特征 2015.6.24
            feature_fmj{t}{j,mm} = [ diff_intensity_sum; diff_size_sum; distance;...
                diff_son; father_fe ];
        end
    end
end
toc

%% 将特征归一化到【-1，1】区间，方便后续计算(2015.9.29)
tic
disp('进行特征归一化...');
if strcmp(dataset, 'training')
    [ feature_fij feature_fij_p minmax.fij.min minmax.fij.max ] = CX_mapminmax( feature_fij );
    [ feature_fid feature_fid_p minmax.fid.min minmax.fid.max ] = CX_mapminmax( feature_fid );
    [ feature_fit feature_fit_p minmax.fit.min minmax.fit.max ] = CX_mapminmax( feature_fit );
    [ feature_fiv feature_fiv_p minmax.fiv.min minmax.fiv.max ] = CX_mapminmax( feature_fiv );
    [ feature_fmj feature_fmj_p minmax.fmj.min minmax.fmj.max ] = CX_mapminmax( feature_fmj );
    [ feature_fsj feature_fsj_p minmax.fsj.min minmax.fsj.max ] = CX_mapminmax( feature_fsj );
    save([ trackpath, '\结构化学习\minmax.mat'], 'minmax');
else
    [ ~, traintrackpath ] = getpath( 'training' ); % 测试需要从训练中载入归一参数
    load([ traintrackpath, '\结构化学习\minmax.mat'], 'minmax');
    [ feature_fij feature_fij_p ] = CX_mapminmax( feature_fij, minmax.fij );
    [ feature_fid feature_fid_p ] = CX_mapminmax( feature_fid, minmax.fid );
    [ feature_fit feature_fit_p ] = CX_mapminmax( feature_fit, minmax.fit );
    [ feature_fiv feature_fiv_p ] = CX_mapminmax( feature_fiv, minmax.fiv );
    [ feature_fmj feature_fmj_p ] = CX_mapminmax( feature_fmj, minmax.fmj );
    [ feature_fsj feature_fsj_p ] = CX_mapminmax( feature_fsj, minmax.fsj );
end
toc

%% 保存特征数据
disp('  保存数据');
% 原始特征只有在计算 SVM 的时候有用
save(Feature_New_addr, 'feature_fid','feature_fij','feature_fit','feature_fiv','feature_fmj','feature_fsj');
% 增广特征另存一下
save(Feature_Plus_New_addr, 'feature_fid_p','feature_fij_p','feature_fit_p','feature_fiv_p','feature_fmj_p','feature_fsj_p');

% 总维数:    ij   id  sj   it    iv  mj
d_feature = 12 + 24 + 11 + 11 + 23 + 23;

