% function CXSL_Combine_feature_all_New( str )
% =========================================================================
% 2015.6.16
% 这个函数与 CXSL_Combine_feature 的区别在于：需要计算全体流程变量对应的特征
% CXSL_Combine_feature 对概率进行了过滤，只需计算较少地区的特征，但在后续的学习中可能出问题
%
% input: Ellipse  是椭圆数组
%        n  是各帧中的椭圆数目     
%
% output: 将各个事件的特征都保存在 \01_2-16_track\结构化学习\Feature.mat 中
%         'feature_ffd','feature_fid','feature_fij','feature_fit','feature_fiv','feature_fmj','feature_fsj'
%
% =========================================================================
clear;close all
if 1
    dataset = 'competition';
else
    dataset = 'training';
end
[ segpath trackpath ] = getpath( dataset );

pre_data_addr = [ trackpath, '\Pair\Pre_data_New.mat'];

last = max(strfind(trackpath, '\'));
rawpic_addr = trackpath(1:last+2);  % 需要给出原始图片的地址
% 原始维度特征
mkdir([trackpath, '\结构化学习']);
Feature_New_addr = [ trackpath, '\结构化学习\Feature_New.mat'];
% 增广特征
Feature_Plus_New_addr = [ trackpath, '\结构化学习\Feature_Plus_New.mat']; 

%% 需要调用 CX_ILP_Pair_Pre 处理好的椭圆数据
load( pre_data_addr );
frame = numel(Ellipse);
% frame = 60;

%% 计算单个椭圆的特征，加入到自身的struct中
disp('计算椭圆自身特征...');
tic;
Ellipse = CXSL_Calculate_Ellipse_feature( Ellipse, rawpic_addr );
toc;

%% 下面开始计算前后2帧椭圆特征差异
disp('计算各细胞事件特征...');

%% ================== move ================== 迁移特征全是特征差的形式
%
% 特征1-4   ：position，size，eccentric，alpha差异
% 特征5     ：16bin 灰度直方图差异 使用了EMD距离 （获得了补充材料后进行的更新 2015.6.8）
% dist= emd_hat_gd_metric_mex([0.5;0.3;0.2], [0.1;0.1;0.8] ,ones(3,3) ,0 ,1)
% 特征6-8   ：灰度和，灰度均值，灰度标准差的差异
%###########################################  时间花费 32秒 
tic;
feature_fij = cell(frame-1,1);
feature_fij_p = cell(frame-1,1);
for t=1:frame-1
    for j=1:n(t)
        
        e_j = Ellipse{t}{j};
        fj_1_4 = struct2array( e_j.feature.geo ); % j单点特征的临时变量
        for mm=1:4
            % 映射寻找 k
            k = candidate_fij{t}(j,mm);
            
            % 候选椭圆e_k
            e_k = Ellipse{t+1}{k};
            fk_1_4 = struct2array( e_k.feature.geo ); % k单点特征的临时变量
            % 特征1-4的差异
            diff_position = norm( fj_1_4(1:2)-fk_1_4(1:2) ); % 位置差异使用欧氏距离
            diff_1_4 = [ diff_position;( fj_1_4(3:end) - fk_1_4(3:end) )']; % size等差异之间相减
            % 计算灰度直方图的EMD
            hist_j = e_j.feature.intensity.hist;
            hist_k = e_k.feature.intensity.hist;
             % --- 调用外部函数 emd_hat_gd_metric_mex 计算 EM 距离
            diff_5 = emd_hat_gd_metric_mex(hist_j, hist_k, ones(numel(hist_j), numel(hist_j)) ,0 ,1);
            % 计算特征6-8的差异
            fj_6_8 = [ e_j.feature.intensity.sum, e_j.feature.intensity.mean, e_j.feature.intensity.devia ];
            fk_6_8 = [ e_k.feature.intensity.sum, e_k.feature.intensity.mean, e_k.feature.intensity.devia ];
            diff_6_8 = (fj_6_8 - fk_6_8)';
            % 合成move特征，最后补上一个1，作为增广特征 2015.6.24（后来移到mapmimax归一化中进行增广的处理）
            feature_fij{t}{j,mm} = [ diff_1_4; diff_5; diff_6_8 ];

        end      
%         others = setdiff( 1:n(t+1),candidate_k );
%         for ind=1:numel(others)
%             feature_fij{t}{j,others(ind)} = Inf;
%         end   
    end

end
toc;

%% ==================  divide ================== 
%
% 特征1-2  ：灰度均值差异，角度信息
% 特征3-7  ：2子细胞ec差异、size差异、灰度均值的差异、shpae_compactness（2细胞面积之和/凸包多边形面积）
% 特征8-9   ：父ec（母细胞通常是狭长的）、父灰度均值
%########################################### 
tic;
feature_fid = cell(frame-1,1);
feature_fid_p = cell(frame-1,1);
for t=1:frame-1
    for j=1:n(t)
        for mm=1:6
            
            sons = candidate_k_next{t}{j,mm};
            son1 = Ellipse{t+1}{sons(1)};
            son2 = Ellipse{t+1}{sons(2)};
            father = Ellipse{t}{j};
            
            % 灰度和差异
            diff_intensity_sum = father.feature.intensity.mean - son1.feature.intensity.sum - son2.feature.intensity.sum;
            % 角度，使用向量内积来求角度
            v1 = son1.feature.geo.position - father.feature.geo.position;
            v2 = son2.feature.geo.position - father.feature.geo.position;
            angle_patten = acosd(dot(v1,v2)/( norm(v1)*norm(v2) ));
            
            % 2个子细胞的ec差异
            diff_sons_ec = abs( son1.feature.geo.eccentric - son2.feature.geo.eccentric );
            % 2个子细胞的size差异比上较小的那个size
            diff_sons_size = abs( son1.feature.geo.size - son2.feature.geo.size )/...
                min(son1.feature.geo.size, son2.feature.geo.size);
            % 2个子细胞的灰度均值差异
            diff_sons_intensity = abs( son1.feature.intensity.mean - son2.feature.intensity.mean );
            % 细胞的shape compactness指标
            shape_compact = pi*(son1.a*son1.b + son2.a*son2.b) /CX_Convhull(son1, son2);
            
            % 父ec
            father_ec = father.feature.geo.eccentric;
            % 父灰度均值
            father_intensity = father.feature.intensity.mean;
            % 合成分裂特征，最后补上一个1，作为增广特征 2015.6.24
            feature_fid{t}{j,mm} = [ diff_intensity_sum, angle_patten, diff_sons_ec, diff_sons_size,...
                diff_sons_intensity, shape_compact, father_ec, father_intensity ]';
        end
    end
end
toc;   

%% ==================  appear ================== 
%
% 特征1-2  ：size、dist2border
% 特征3-5  ：灰度和，灰度均值，灰度标准差
%
%########################################### 
tic
feature_fsj = cell(frame-1,1);
feature_fsj_p = cell(frame-1,1);
for t=2:frame
    for j=1:n(t)
        e_j = Ellipse{t}{j};
        s_d = [ e_j.feature.geo.size, e_j.feature.dist2border ];
        feature_intensity = [ e_j.feature.intensity.sum, e_j.feature.intensity.mean, e_j.feature.intensity.devia ];
        % 最后补上一个1，作为增广特征 2015.6.24
        feature_fsj{t}{j,1} = [ s_d, feature_intensity ]';
    end
end
toc

%% ==================  disappear ================== 
%
% 特征1-2  ：size、dist2border
% 特征3-5  ：灰度和，灰度均值，灰度标准差
%
%########################################### 
tic
feature_fit = cell(frame-1,1);
feature_fit_p = cell(frame-1,1);
for t=1:frame-1
    for j=1:n(t)
        e_j = Ellipse{t}{j};
        s_d = [ e_j.feature.geo.size, e_j.feature.dist2border ];
        feature_intensity = [ e_j.feature.intensity.sum, e_j.feature.intensity.mean, e_j.feature.intensity.devia ];
         % 最后补上一个1，作为增广特征 2015.6.24
        feature_fit{t}{j,1} = [ s_d, feature_intensity ]';
    end
end
toc

%% ==================  split ================== 
%
% 特征1-2  ：灰度和差异、size和差异
% 特征3-4  ：子细胞shape_compactness、子细胞size差异
% 特征5    ：源细胞的拟合度指标hd（删去）2015.9.29
%########################################### 
tic
feature_fiv = cell(frame-1,1);
feature_fiv_p = cell(frame-1,1);
for t=1:frame-1
    for j=1:n(t)
        for mm=1:6
            
            couples = candidate_k_next{t}{j,mm};
            brother1 = Ellipse{t+1}{couples(1)};
            brother2 = Ellipse{t+1}{couples(2)};
            father = Ellipse{t}{j};
            
            % 灰度和差异
            diff_intensity_sum = father.feature.intensity.sum - brother1.feature.intensity.sum - brother2.feature.intensity.sum;
            % size和差异
            diff_size_sum = father.feature.geo.size - brother1.feature.geo.size - brother2.feature.geo.size;
            
            % 细胞的shape compactness指标
            shape_compact = pi*(brother1.a*brother1.b + brother2.a*brother2.b) / CX_Convhull(brother1, brother2);
            % 子细胞size差异
            diff_sons_size = abs( brother1.feature.geo.size - brother2.feature.geo.size );
            
            % 源细胞的拟合度指标hd
%             source_fitness = father.hd;
            % 最后补上一个1，作为增广特征 2015.6.24
            feature_fiv{t}{j,mm} = [ diff_intensity_sum, diff_size_sum, shape_compact, diff_sons_size ]';
        end
    end
end
toc

%% ==================  merge ================== 
%
% 特征1-2  ：灰度和差异、size和差异
% 特征3-4  ：子细胞shape_compactness、源细胞size差异
% 特征5    ：大细胞的拟合度指标hd（删去）2015.9.29
%########################################### 
tic
feature_fmj = cell(frame-1,1);
feature_fmj_p = cell(frame-1,1);
for t=2:frame
    for j=1:n(t)
        for mm=1:6

            sources = candidate_k_last{t}{j,mm};
            source1 = Ellipse{t-1}{sources(1)};
            source2 = Ellipse{t-1}{sources(2)};
            father = Ellipse{t}{j};
            
            % 灰度和差异
            diff_intensity_sum = father.feature.intensity.sum- source1.feature.intensity.sum - source2.feature.intensity.sum;
            % size和差异
            diff_size_sum = father.feature.geo.size - source1.feature.geo.size - source2.feature.geo.size;
            
            % 细胞的shape compactness指标
            shape_compact = pi*(source1.a*source1.b + source2.a*source2.b) / CX_Convhull(source1, source2);
            % 源细胞size差异
            diff_sons_size = abs( source1.feature.geo.size - source2.feature.geo.size );
            
            % 大细胞的拟合度指标hd
%             source_fitness = father.hd;
            % 最后补上一个1，作为增广特征 2015.6.24
            feature_fmj{t}{j,mm} = [ diff_intensity_sum, diff_size_sum, shape_compact, diff_sons_size ]';
        end
    end
end
toc

%% ==================  false detetion ================== 
%
% 特征1     ：size
% 特征2-4   ：灰度和，灰度均值，灰度标准差, dist2border
%###########################################
feature_ffd = cell(frame-1,1);
feature_ffd_p = cell(frame-1,1);
for t=1:frame-1
    for j=1:n(t)
        p_s = Ellipse{t}{j}.feature.geo.size;
        feature_intensity = [ Ellipse{t}{j}.feature.intensity.sum, Ellipse{t}{j}.feature.intensity.mean, Ellipse{t}{j}.feature.intensity.devia ];
        % 最后补上一个1，作为增广特征 2015.6.24
        feature_ffd{t}{j,1} = [ p_s, feature_intensity, Ellipse{t}{j}.feature.dist2border ]';
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

     
        
 
%% 这个函数用于计算2个椭圆的最小凸包多边形的面积
% 主要利用内置函数 convhull 进行计算
if 0 
% function Area = CX_Convhull(e1, e2)
     
% 先计算e1的圆周点，取361即可
alpha1 = e1.alpha;
a = e1.a;
b = e1.b;
x0 = e1.x0;
y0 = e1.y0;

c=cosd(alpha1);
s=sind(alpha1);
polar_angle=linspace(0,360,50);
xq= a*cosd(polar_angle);
yq= b*sind(polar_angle);
xn1=xq*c-yq*s+x0;
yn1=xq*s+yq*c+y0;

% 再计算e2的圆周点，取20即可
alpha1 = e2.alpha;
a = e2.a;
b = e2.b;
x0 = e2.x0;
y0 = e2.y0;

c=cosd(alpha1);
s=sind(alpha1);
polar_angle=linspace(0,360,50);
xq= a*cosd(polar_angle);
yq= b*sind(polar_angle);
xn2=xq*c-yq*s+x0;
yn2=xq*s+yq*c+y0;

xx = [xn1 xn2];
yy = [yn1 yn2];

[~, Area] = convhull(xx, yy);       
        
end

