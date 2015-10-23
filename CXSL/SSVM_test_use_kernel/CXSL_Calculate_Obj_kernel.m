function object_function = CXSL_Calculate_Obj_kernel( dataset, s_frame, e_frame,...
    A_best, alpha_best, kernel_type, cmd,...
    fij, fit, fid, fiv, fmj, fsj )

% ======================================================================= %
%
% 这个函数仅用于debug，其函数原型是作为 CXSL_ILP_Using_Best_W 内部调用
% 功能是实现另一种计算目标函数的方法，即 <w,feature>.*z
%
% ======================================================================= %

% ----------------------------------------------------------------------- %
% 载入特征
[ ~, trackpath ] = getpath( dataset );
load([ trackpath, '\结构化学习\Feature_New.mat']);
% ----------------------------------------------------------------------- %

obj_ij = 0;
obj_it = 0;
obj_id = 0;
obj_iv = 0;
obj_mj = 0;
obj_sj = 0;

for tt=s_frame:e_frame-1
    %% ------------------------- fij -------------------------- %
    ev = 1;
    n_A = size(A_best{ev},2);
    n_F = numel(feature_fij{tt});    
    K = zeros(n_A, n_F);
    for i=1:n_A
        for j=1:n_F         
            K(i,j) = svm_kernel(A_best{ev}(:,i), feature_fij{tt}{j}, kernel_type{ev}, cmd{ev});
        end
    end
    % 当前帧的目标函数为(1,n_A) * (n_A,n_F) * (n_F,1) 
    obj_ij = obj_ij + alpha_best{ev}'*K*reshape(fij{tt},[],1);
    % ------------------------------------------------------- %
    
    %% ------------------------- fit -------------------------- %
    ev = 2;
    n_A = size(A_best{ev},2);
    n_F = numel(feature_fit{tt});    
    K = zeros(n_A, n_F);
    for i=1:n_A
        for j=1:n_F         
            K(i,j) = svm_kernel(A_best{ev}(:,i), feature_fit{tt}{j}, kernel_type{ev}, cmd{ev});
        end
    end
    % 当前帧的目标函数为(1,n_A) * (n_A,n_F) * (n_F,1) 
    obj_it = obj_it + alpha_best{ev}'*K*reshape(fit{tt},[],1);
    % ------------------------------------------------------- %
    
    %% ------------------------- fid -------------------------- %
    ev = 3;
    n_A = size(A_best{ev},2);
    n_F = numel(feature_fid{tt});    
    K = zeros(n_A, n_F);
    for i=1:n_A
        for j=1:n_F         
            K(i,j) = svm_kernel(A_best{ev}(:,i), feature_fid{tt}{j}, kernel_type{ev}, cmd{ev});
        end
    end
    % 当前帧的目标函数为(1,n_A) * (n_A,n_F) * (n_F,1) 
    obj_id = obj_id + alpha_best{ev}'*K*reshape(fid{tt},[],1);
    % ------------------------------------------------------- %
    
    %% ------------------------- fiv -------------------------- %
    ev = 4;
    n_A = size(A_best{ev},2);
    n_F = numel(feature_fiv{tt});    
    K = zeros(n_A, n_F);
    for i=1:n_A
        for j=1:n_F         
            K(i,j) = svm_kernel(A_best{ev}(:,i), feature_fiv{tt}{j}, kernel_type{ev}, cmd{ev});
        end
    end
    % 当前帧的目标函数为(1,n_A) * (n_A,n_F) * (n_F,1) 
    obj_iv = obj_iv + alpha_best{ev}'*K*reshape(fiv{tt},[],1);
    % ------------------------------------------------------- %
    
end
for tt=s_frame+1:e_frame
    %% ------------------------- fmj -------------------------- %
    ev = 5;
    n_A = size(A_best{ev},2);
    n_F = numel(feature_fmj{tt});    
    K = zeros(n_A, n_F);
    for i=1:n_A
        for j=1:n_F         
            K(i,j) = svm_kernel(A_best{ev}(:,i), feature_fmj{tt}{j}, kernel_type{ev}, cmd{ev});
        end
    end
    % 当前帧的目标函数为(1,n_A) * (n_A,n_F) * (n_F,1) 
    obj_iv = obj_iv + alpha_best{ev}'*K*reshape(fmj{tt},[],1);
    % ------------------------------------------------------- %
    
    %% ------------------------- fsj -------------------------- %
    ev = 6;
    n_A = size(A_best{ev},2);
    n_F = numel(feature_fsj{tt});    
    K = zeros(n_A, n_F);
    for i=1:n_A
        for j=1:n_F         
            K(i,j) = svm_kernel(A_best{ev}(:,i), feature_fsj{tt}{j}, kernel_type{ev}, cmd{ev});
        end
    end
    % 当前帧的目标函数为(1,n_A) * (n_A,n_F) * (n_F,1) 
    obj_sj = obj_sj + alpha_best{ev}'*K*reshape(fsj{tt},[],1);
    % ------------------------------------------------------- %
    
end
    
object_function = obj_ij + obj_it + obj_id + obj_iv + obj_mj + obj_sj;

end



