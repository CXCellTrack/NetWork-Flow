function [ fai_x_z ] = CXSL_Calculate_phi_x_zstar_New( w, s_frame, e_frame, str)

% ============================= %
% 这个函数用于计算 fai(x,z*) 
% 也可以计算 fai(x,z^) 但 fai(x,z^) 可以由ILP中的结果 value 出来，因此无需重复计算
% ============================= %

if isempty(w) % 若没有输入w则默认为带b的0
    w = zeros(40,1);
end

dataset = 'training';
[ ~, trackpath ] = getpath( dataset );
% 载入特征
if numel(w)==40
    load([ trackpath, '\结构化学习\Feature_Plus_New.mat']);
    feature_fij = feature_fij_p;
    feature_fit = feature_fit_p;
    feature_fid = feature_fid_p;
    feature_fiv = feature_fiv_p;
    feature_fmj = feature_fmj_p;
    feature_fsj = feature_fsj_p;
else
    load([ trackpath, '\结构化学习\Feature_New.mat']);
end
% 载入标准答案GT
load([ trackpath, '\GT\GT_after_hand_tune\GT_Flow_Variables_New.mat']);

fai_fij = 0;
fai_fit = 0;
fai_fid = 0;
fai_fiv = 0;
fai_fmj = 0;
fai_fsj = 0;

if strcmp(str, 'hat') % 这部分已经由ILP的结果中value出，无需重新计算
    %
end
if strcmp(str, 'star')
    % 计算标准答案中的fai(x,z*)
%     tic;
    for t = s_frame:e_frame-1  
        
        fai_fij = fai_fij + cell_dot_mutil( feature_fij{t}, Fij{t} );
            
        fai_fit = fai_fit + cell_dot_mutil( feature_fit{t}, Fit{t} ); %  加不加sum = 都一样
        
        fai_fid = fai_fid + cell_dot_mutil( feature_fid{t}, Fid{t} );  
		
        fai_fiv = fai_fiv + cell_dot_mutil( feature_fiv{t}, Fiv{t} );         
    end
    
    for t = s_frame+1:e_frame
	
        fai_fmj = fai_fmj + cell_dot_mutil( feature_fmj{t}, Fmj{t} );
		
        fai_fsj = fai_fsj + cell_dot_mutil( feature_fsj{t}, Fsj{t} );
    end
%     toc;
end

% 将所有的 fai_f 组合成列向量（增广向量）一共42维
fai_x_z = [ fai_fij; fai_fit; fai_fid; fai_fiv; fai_fmj; fai_fsj ];

end


%% 这个函数与 CXSL_Calculate_Obj 中的函数一样
% 注意：由于数据为double类型，使用for循环反而快

function sum_f = cell_dot_mutil( feature, z )
%
% f 为特征向量cell	n(t)*n(t+1)
% g 为流程变量矩阵
% 输出 sum_f 为该cell内所有特征向量之和
% ===================================

% 1. for 循环式计算 0.035 秒

sum_f = zeros( size(feature{1,1}) );
[h w] = size(feature);
for i=1:h
    for j=1:w
        sum_f = sum_f + feature{i, j}* z(i,j);
    end
end

% 2. arrayfun运算 1.77 秒

% sum_f = 0;
% ss = numel(feature);
% % feature1 = cellfun(@(x)x', feature, 'un',0);
% feature2 = reshape(feature, ss, 1);
% z1 = reshape(z, ss, 1);
% f_z = arrayfun(@(x)feature2{x}*z1(x), 1:ss, 'un',0);
% for x=1:ss
%     sum_f = sum_f + f_z{x};
% end

end
















