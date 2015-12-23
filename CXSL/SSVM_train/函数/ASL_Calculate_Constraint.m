function F = ASL_Calculate_Constraint( dataset, s_frame,e_frame,fij,fit,fid,fiv,fmj,fsj )
% ================================================================== %
%
% 本函数用于计算约束条件 F，active sl方法专用
% 调用了 CXSL_Assign_FlowVar_With_Loss 来预分配流程变量以及计算损失函数
% 最终输出各流程变量、损失函数以及约束条件
%
% ASL的约束与OURS的区别：
%       1、ASL给每个前景分配一个假说，无论是单假说前景还是多假说前景
%       2、ASL的每个假说都必须被解释
%
% ================================================================== %
global Ellipse n conflict_fij conflict_pair_next_xy conflict_pair_last_xy;

[ ~, trackpath ] = getpath( dataset );
% 载入数据，可以选择载入训练集或测试集上的数据

if exist('SuperPixel','var') % 如果用的是超像素，就采用赋值
    disp('  采用的是SuperPixel假说！');
    Ellipse = SuperPixel;
    clear SuperPixel
end

F1 = [];
F2 = [];
F3 = [];
F4 = [];

%######################################## 建立约束条件 #############################################
disp('   构建约束1：进出约束……');
assert(s_frame+1==e_frame); % 只使用于2帧的情况

%% 约束1
%################################################################
% 针对数据集3（Fluo-N2DH-SIM+）的约束（禁止merge和split事件的发生！椭圆假说中，超像素不禁止）
if ~isempty(strfind(trackpath,'SIM+'))
    disp('Attention! merge&split event has been canceled！')
    for t = s_frame+1:e_frame
        F1 = [ F1, sum(fmj{t}(:))==0 ];
    end
    for t = s_frame:e_frame-1
        F1 = [ F1, sum(fiv{t}(:))==0 ];
    end
end

%% 约束2
% 融合之后不允许立刻分裂
for t = s_frame+1:e_frame-1
    for j=1:n(t)
        F2 = [ F2, sum(fmj{t}(j,:)) + sum(fid{t}(j,:)) <= 1 ];
    end
end

%% 约束3：每个单假说前景必须被解释
t = e_frame; % 2-t 帧要求每个前景内至少有一个被解释（入口解释）
for j=1:n(t)
    sum_fid = 0;
    sum_fiv = 0;
    sum_fij = 0;
    % 单独前景必须被解释，即入口和为1
    if Ellipse{t}{j}.num_hypoth == 1
        
        % sum_fij 为所有迁移到 j 的fij之和（入口和），出口和可以用 sum(fij{t}(j,:)) 表示
        for ind=1:size(conflict_fij{t-1}{j}, 1)
            sum_fij = sum_fij + fij{t-1}( conflict_fij{t-1}{j}(ind,1), conflict_fij{t-1}{j}(ind,2) );
        end
        
        % sum_fid 为所有分裂到包含 j 的 pair 的 fid 之和
        for ind=1:numel(conflict_pair_next_xy{t}{j})/2
            sum_fid = sum_fid + fid{t-1}( conflict_pair_next_xy{t}{j}(ind,1), conflict_pair_next_xy{t}{j}(ind,2) );
            sum_fiv = sum_fiv + fiv{t-1}( conflict_pair_next_xy{t}{j}(ind,1), conflict_pair_next_xy{t}{j}(ind,2) );
        end
        F3 = [ F3, sum_fij + fsj{t}(j) + sum(fmj{t}(j,:)) + sum_fid + sum_fiv == 1 ];
    else
        % 多目标前景中至少有一个要被解释
        %
        % 这部分写在j到j+jplus为一个前景处
    end
end

% 同时需要添加1~t-1帧的出口必须被解释
t = s_frame;
for j=1:n(t)
    if Ellipse{t}{j}.num_hypoth == 1
        sum_fmj = 0;
        % sum_fmj 为所有包含 j 的融合 pair 的 fmj 之和
        for ind=1:numel(conflict_pair_last_xy{t}{j})/2
            sum_fmj = sum_fmj + fmj{t+1}( conflict_pair_last_xy{t}{j}(ind,1), conflict_pair_last_xy{t}{j}(ind,2) );
        end
        F3 = [ F3, sum(fij{t}(j,:)) + fit{t}(j) + sum(fid{t}(j,:)) + sum(fiv{t}(j,:)) + sum_fmj == 1 ];
    else
        % 多目标前景中至少有一个要被解释
        %
        % 这部分写在j到j+jplus为一个前景处
    end
end

disp('   构建约束2：多假说约束……');

%% 约束4：每个多假说前景中仅有一个成立

for t = s_frame:e_frame
    j = 1;
    while j<=n(t)    % j遍历
        %% 单前景假说则j自增1，继续寻找下个多假说前景
        if Ellipse{t}{j}.num_hypoth==1
            j = j + 1;
            continue;
        end
        
        % 这一部分求出j到j+jplus都属于同一前景
        jplus = Ellipse{t}{j}.num_hypoth - 1;
        % 下面这段实现相同的功能，但代码复杂
        %         jplus = 1;
        %         % j+jplus不能超过总数
        %         while j+jplus<=numel(Ellipse{t}) && Ellipse{t}{j+jplus}.ind_region == Ellipse{t}{j}.ind_region
        %             jplus = jplus + 1;
        %         end
        %         jplus = jplus - 1;
        %##################### j到j+jplus为一个前景 ####################

        %% 此处添加一个多目标前景必须要有一个入口的条件
        %##########################################
        % 第一帧的多目标前景出口为1
        if t==s_frame
            % 累加起来
            sum_out = 0;
            for thisj=j:j+jplus
                sum_fmj = 0;
                % 求出前景内所有细胞merge出口的和
                for ind=1:numel(conflict_pair_last_xy{t}{thisj})/2
                    sum_fmj = sum_fmj + fmj{t+1}( conflict_pair_last_xy{t}{thisj}(ind,1), conflict_pair_last_xy{t}{thisj}(ind,2) );
                end
                % 求出前景内所有细胞其他 4 种出口的和
                sum_fij = sum(fij{t}(thisj,:));
                sum_fit = fit{t}(thisj);
                sum_fid = sum(fid{t}(thisj,:));
                sum_fiv = sum(fiv{t}(thisj,:));
                
                sum_out = sum_out + sum_fmj + sum_fij + sum_fit + sum_fid + sum_fiv;
                
            end
            F4 = [ F4,  sum_out==1 ];
        else
            % 第二帧的多目标前景入口为1
            sum_in = 0;
            for thisj=j:j+jplus
                
                sum_fij = 0;
                sum_fid = 0;
                sum_fiv = 0;
                % 求出前景内所有细胞 split 和 divide 入口的和
                for ind=1:numel(conflict_pair_next_xy{t}{thisj})/2
                    sum_fid = sum_fid + fid{t-1}( conflict_pair_next_xy{t}{thisj}(ind,1), conflict_pair_next_xy{t}{thisj}(ind,2) );
                    sum_fiv = sum_fiv + fiv{t-1}( conflict_pair_next_xy{t}{thisj}(ind,1), conflict_pair_next_xy{t}{thisj}(ind,2) );
                end
                % 求出前景内所有细胞其他 3 种入口的和

                % sum_fij 为所有迁移到 j 的fij之和（入口和），出口和可以用 sum(fij{t}(j,:)) 表示
                for ind=1:size(conflict_fij{t-1}{thisj}, 1)
                    sum_fij = sum_fij + fij{t-1}( conflict_fij{t-1}{thisj}(ind,1), conflict_fij{t-1}{thisj}(ind,2) );
                end

                sum_fmj = sum(fmj{t}(thisj,:));
                sum_fsj = fsj{t}(thisj);
                
                sum_in = sum_in + sum_fmj + sum_fij + sum_fsj + sum_fid + sum_fiv;
                
            end
            F4 = [ F4,  sum_in==1 ];
        end
        %##########################################
        %
        % 以上一段是在一个多目标前景内设置入口至少为1的约束
        % 即保证每个多目标前景都必须被解释，不允许虚景的存在
        %
        %##############################################################
       j = j + jplus + 1; % j直接走出这个前景 
       
    end  %## end while
end

%% 导出最终的目标函数与总约束

%#################################
% 各个约束条件的解释如下：
%######### 必要约束 4 个 ##########

% F1：针对数据集3（Fluo-N2DH-SIM+）的约束（禁止merge和split事件的发生！
% F2：融合之后不允许立刻分裂
% F3：每个单假说前景必须被解释
% F4：每个多假说前景中仅有一个成立


F = [ F1, F2, F3, F4 ];

end



