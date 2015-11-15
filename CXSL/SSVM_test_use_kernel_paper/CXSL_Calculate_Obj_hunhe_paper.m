function object_function = CXSL_Calculate_Obj_hunhe_paper( dataset, s_frame, e_frame,...
    w_best, delta_y_best, feature_best, alpha_best, kernel_type, cmd, lambda, islinear, ...
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
    if islinear(ev)
        mat_fij = cellfun(@(x)dot(w_best{ev},x), feature_fij{tt});
        obj_ij = obj_ij + sum(sum( mat_fij.*fij{tt} ));
    end
    % ------------------------------------------------------- %
    
    %% ------------------------- fit -------------------------- %
    ev = 2;
    if islinear(ev)
        mat_fit = cellfun(@(x)dot(w_best{ev},x), feature_fit{tt});
        obj_it = obj_it + sum(sum( mat_fit.*fit{tt} ));
    end
    % ------------------------------------------------------- %
    
    %% ------------------------- fid -------------------------- %
    ev = 3;
    if islinear(ev)
        mat_fid = cellfun(@(x)dot(w_best{ev},x), feature_fid{tt});
        obj_id = obj_id + sum(sum( mat_fid.*fid{tt} ));
    end
    % ------------------------------------------------------- %
    
    %% ------------------------- fiv -------------------------- %
    ev = 4;
    if islinear(ev)
        mat_fiv = cellfun(@(x)dot(w_best{ev},x), feature_fiv{tt});
        obj_iv = obj_iv + sum(sum( mat_fiv.*fiv{tt} ));
    end    
    % ------------------------------------------------------- %
    
end
for tt=s_frame+1:e_frame
    %% ------------------------- fmj -------------------------- %
    ev = 5;
    if islinear(ev)
        mat_fmj = cellfun(@(x)dot(w_best{ev},x), feature_fmj{tt});
        obj_mj = obj_mj + sum(sum( mat_fmj.*fmj{tt} ));
    end    
    % ------------------------------------------------------- %
    
    %% ------------------------- fsj -------------------------- %
    ev = 6;
    if islinear(ev)
        mat_fsj = cellfun(@(x)dot(w_best{ev},x), feature_fsj{tt});
        obj_sj = obj_sj + sum(sum( mat_fsj.*fsj{tt} ));
    end    
    % ------------------------------------------------------- %
    
end
    
for ev=1:6
    if islinear(ev)
        continue;
    end
    %% 考察分裂事件
%     tic
    if ev==3 
        N = size(alpha_best,1);
        for ii=1:N
            for nsv=1:numel(alpha_best{ii,ev})
                alpha1 = alpha_best{ii,ev}(nsv); % 一个alpha
                if alpha1==0
                    continue;
                end
                
                disp(['样本',num2str(ii),',支持向量',num2str(nsv)]);
                delta_y1 = delta_y_best{ii,ev}(nsv,:); % 对应的一个样本delta_y
                feature1 = feature_best{ii,ev};
                
                frameK = 0; % 累加所有帧
                for tt1=1:numel(delta_y1) % 按帧走
                    for tt2=1:numel(fid) % 按帧走
%                         disp(['样本',num2str(ii),',支持向量',num2str(nsv),',计算前者第',num2str(tt1),'帧与后者第',num2str(tt2),'帧...']);
                        y1 = delta_y1{tt1};
                        y2 = fid{tt2};
                        
                        K_mat = zeros(numel(y1), numel(y2));
                        for kk=1:numel(y1)
                            if y1(kk)==0
                                continue;
                            end
                            for jj=1:numel(y2)
                                K_mat(kk,jj) = svm_kernel(feature1{tt1}{kk}, feature_fid{tt2}{jj}, kernel_type{ev}, cmd{ev});
                            end
                        end
                        
                        ytmp1 = reshape(y1,[],1)';
                        ytmp2 = reshape(y2,[],1);
                        frameK = frameK + ytmp1*K_mat*ytmp2;
                    end
                end
                
                obj_id = obj_id + alpha1*frameK;
            end
        end
    end
    
    obj_id = obj_id/(lambda(ev)*N); % lambda得加上
%     toc
    
    %% 别的事件先不考察了
    
end
              
                

object_function = obj_ij + obj_it + obj_id + obj_iv + obj_mj + obj_sj;




























end



