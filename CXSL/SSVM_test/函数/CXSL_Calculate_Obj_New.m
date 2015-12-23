function object_function = CXSL_Calculate_Obj_New( dataset, w_best, s_frame, e_frame, fij, fit, fid, fiv, fmj, fsj )

% ======================================================================= %
%
% 这个函数仅用于debug，其函数原型是作为 CXSL_ILP_Using_Best_W 内部调用
% 功能是实现另一种计算目标函数的方法，即 <w,feature>.*z
%
% ======================================================================= %

% ----------------------------------------------------------------------- %
[ ~, trackpath ] = getpath( dataset );
% 进来的 w_best 为列向量，此处先转换为行向量再进行相乘
% 载入特征
if numel(w_best)>=40
    load([ trackpath, '\结构化学习\Feature_Plus_New.mat']);
else
    load([ trackpath, '\结构化学习\Feature_New.mat']);
end

if exist('feature_fij_p','var')
    % 如为增广的特征则使用增广的w 
    d_f = [ 13, 12, 25, 24, 24, 12 ]; % 各事件特征的维数(要按顺序)
    d_f = [ 9, 6, 9, 5, 5, 6 ];
    tmp = cumsum(d_f);
    wij = w_best( 1:tmp(1) )'; 
    wit = w_best( tmp(1)+1:tmp(2) )';
    wid = w_best( tmp(2)+1:tmp(3) )';
    wiv = w_best( tmp(3)+1:tmp(4) )';
    wmj = w_best( tmp(4)+1:tmp(5) )';
    wsj = w_best( tmp(5)+1:tmp(6) )'; 

    feature_fij = feature_fij_p;
    feature_fit = feature_fit_p;
    feature_fid = feature_fid_p;
    feature_fiv = feature_fiv_p;
    feature_fmj = feature_fmj_p;
    feature_fsj = feature_fsj_p;
else
    % 不含增广的w
    wij = w_best( 1:8 )'; 
    wit = w_best( 9:13 )';
    wid = w_best( 14:21 )';
    wiv = w_best( 22:25 )';
    wmj = w_best( 26:29 )';
    wsj = w_best( 30:34 )';
   
end
    
% 组合目标函数
% ----------------------------------------------------------------------- %

obj_ij = 0;
obj_it = 0;
obj_id = 0;
obj_iv = 0;
obj_mj = 0;
obj_sj = 0;

for tt=s_frame:e_frame-1
    % <w,feature>刚好得到一个矩阵，与fij点乘即可，速度加快
    mat_fij = cellfun(@(x)wij*x, feature_fij{tt});
    obj_ij = obj_ij + sum(sum( mat_fij.*fij{tt} ));
    
    mat_fit = cellfun(@(x)wit*x, feature_fit{tt});
    obj_it = obj_it + sum(sum( mat_fit.*fit{tt} ));
    
    mat_fid = cellfun(@(x)wid*x, feature_fid{tt});
    % 这里有些特征为nan，看来是椭圆假说出了问题！先用下面这句弥补下
    mat_fid(isnan(mat_fid)) = 0;
    obj_id = obj_id + sum(sum( mat_fid.*fid{tt} ));
    
    mat_fiv = cellfun(@(x)wiv*x, feature_fiv{tt});
    obj_iv = obj_iv + sum(sum( mat_fiv.*fiv{tt} ));
end
for tt=s_frame+1:e_frame
    
    mat_fmj = cellfun(@(x)wmj*x, feature_fmj{tt});
    obj_mj = obj_mj + sum(sum( mat_fmj.*fmj{tt} )) ;
    
    mat_fsj = cellfun(@(x)wsj*x, feature_fsj{tt});
    obj_sj = obj_sj + sum(sum( mat_fsj.*fsj{tt} ));
    
end
    
object_function = obj_ij + obj_it + obj_id + obj_iv + obj_mj + obj_sj;

end



