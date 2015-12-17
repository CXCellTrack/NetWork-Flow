function [ Ffull Fbase ] = CXSL_Calculate_Constraint_New_Conflict( dataset, use_op_cons, s_frame,e_frame,fij,fit,fid,fiv,fmj,fsj )
% ================================================================== %
%
% 本函数用于计算约束条件 F
% 调用了 CXSL_Assign_FlowVar_With_Loss 来预分配流程变量以及计算损失函数
% 最终输出各流程变量、损失函数以及约束条件
%
% ================================================================== %

[ ~, trackpath ] = getpath( dataset );
% 载入数据，可以选择载入训练集或测试集上的数据
Pre_data_addr = [ trackpath, '\Pair\Pre_data_New.mat' ];  
load( Pre_data_addr);
if exist('SuperPixel','var') % 如果用的是超像素，就采用赋值
    disp('  采用的是SuperPixel假说！');
    Ellipse = SuperPixel;
    clear SuperPixel
end
    
F1 = [];
F2 = [];
Fop1 = [];
Fop2 = [];
Fop3 = [];
Fop4 = [];
Fop5 = [];
F3 = [];
F4 = [];
Fop6 = [];
Fop7 = [];

%######################################## 建立约束条件 #############################################
disp('   构建约束1：进出守恒……');
tic
%% 约束条件1： 进出守恒

if e_frame - s_frame ~= 1 % 必须在帧数不为2时才有此条件
    disp('      1)一致性约束……');
    %#############  中间帧的入口出口守恒  ###########
    for t = s_frame+1:e_frame-1
        disp(['      正在计算第',num2str(t),'帧的进出守恒...']);
        for j=1:n(t)
            sum_fij = 0;
            sum_fid = 0;
            sum_fiv = 0;
            sum_fmj = 0;

            % sum_fij 为所有迁移到 j 的fij之和（入口和），出口和可以用 sum(fij{t}(j,:)) 表示 
            for ind=1:size(conflict_fij{t-1}{j}, 1)
                sum_fij = sum_fij + fij{t-1}( conflict_fij{t-1}{j}(ind,1), conflict_fij{t-1}{j}(ind,2) );
            end

            % sum_fid 为所有分裂到包含 j 的 pair 的 fid 之和
            for ind=1:numel(conflict_pair_next_xy{t}{j})/2
                sum_fid = sum_fid + fid{t-1}( conflict_pair_next_xy{t}{j}(ind,1), conflict_pair_next_xy{t}{j}(ind,2) );
                sum_fiv = sum_fiv + fiv{t-1}( conflict_pair_next_xy{t}{j}(ind,1), conflict_pair_next_xy{t}{j}(ind,2) );
            end
            % sum_fmj 为所有包含 j 的融合 pair 的 fmj 之和
            for ind=1:numel(conflict_pair_last_xy{t}{j})/2
                sum_fmj = sum_fmj + fmj{t+1}( conflict_pair_last_xy{t}{j}(ind,1), conflict_pair_last_xy{t}{j}(ind,2) );
            end
            F1 = [ F1, sum_fij + fsj{t}(j) + sum(fmj{t}(j,:)) + sum_fid + sum_fiv == sum(fij{t}(j,:)) + fit{t}(j) + sum(fid{t}(j,:)) + sum(fiv{t}(j,:)) + sum_fmj ];
        end
    end
    
end
% end if

%#############  2到t帧入口唯一  ###########
disp('      2)入口约束……');
for t = s_frame+1:e_frame
    for j=1:n(t)
        sum_fid = 0;
        sum_fiv = 0;
        sum_fij = 0;
        % sum_fij 为所有迁移到 j 的fij之和（入口和），出口和可以用 sum(fij{t}(j,:)) 表示 
        for ind=1:size(conflict_fij{t-1}{j}, 1)
            sum_fij = sum_fij + fij{t-1}( conflict_fij{t-1}{j}(ind,1), conflict_fij{t-1}{j}(ind,2) );
        end
        
        % sum_fid 为所有分裂到包含 j 的 pair 的 fid 之和
        for ind=1:numel(conflict_pair_next_xy{t}{j})/2
            sum_fid = sum_fid + fid{t-1}( conflict_pair_next_xy{t}{j}(ind,1), conflict_pair_next_xy{t}{j}(ind,2) );
            sum_fiv = sum_fiv + fiv{t-1}( conflict_pair_next_xy{t}{j}(ind,1), conflict_pair_next_xy{t}{j}(ind,2) );
        end
        % 入口最多有一条   迁移来           新出现      merge融合体      子细胞
        F2 = [ F2, sum_fij + fsj{t}(j) + sum(fmj{t}(j,:)) + sum_fid + sum_fiv <= 1 ];
    end
end

%### 1到t-1帧出口唯一（在一致性约束存在时，可简化为：第一帧出口唯一）  ###%
disp('      3)出口约束……');
for t=s_frame%:e_frame-1
    for j=1:n(t)
        sum_fmj = 0;
        % sum_fmj 为所有包含 j 的融合 pair 的 fmj 之和
        for ind=1:numel(conflict_pair_last_xy{t}{j})/2
            sum_fmj = sum_fmj + fmj{t+1}( conflict_pair_last_xy{t}{j}(ind,1), conflict_pair_last_xy{t}{j}(ind,2) );
        end
        % 出口最多有一条   迁移出去        消失        母细胞        下一帧split        下一帧融合
        F2 = [ F2, sum(fij{t}(j,:)) + fit{t}(j) + sum(fid{t}(j,:)) + sum(fiv{t}(j,:)) + sum_fmj <= 1 ];
    end
end

toc;
disp('   构建约束2：可选约束……');
disp(use_op_cons);
tic

%% 可选约束1
%################################################################
% 针对数据集3（Fluo-N2DH-SIM+）的约束（禁止merge和split事件的发生！椭圆假说中，超像素不禁止）
% if ~isempty(strfind(trackpath,'SIM+'))
%     disp('Attention! merge&split event has been canceled！')
%     for t = s_frame+1:e_frame
%         Fop1 = [ Fop1, sum(fmj{t}(:))<=0 ];
%     end
%     for t = s_frame:e_frame-1
%         Fop1 = [ Fop1, sum(fiv{t}(:))<=0 ];
%     end
% end

if any(use_op_cons==1)
    %################## 可选择的约束一 #################（与真实情况有点矛盾）
    % 分裂出去的不能仍在同一个前景中
    for t = s_frame:e_frame-1
        for j=1:n(t)
            tmp_conflict = 0;
            for mm=1:6
                sons = candidate_k_next{t}{j,mm};
                if Ellipse{t+1}{sons(1)}.ind_region == Ellipse{t+1}{sons(2)}.ind_region % 分裂出去的不能仍在同一个前景中
                    tmp_conflict = tmp_conflict + fid{t}(j,mm);
                end
            end
            % 如果邻域内都是单独前景，tmp_conflict可能为double 0，因此要判断
            if ~isa(tmp_conflict, 'double')
                Fop1 = [ Fop1, tmp_conflict <= 0];
            end
        end
    end
end

%% 可选约束2
if any(use_op_cons==2)
    %################## 可选择的约束二 #################（与真实情况有点矛盾）
    % 不允许不是同一个前景的细胞发生融合 在pair形式下实现这个很容易，只需要比较融合pair中的2细胞是否在一个前景内就可以了
    for t = s_frame+1:e_frame
        for j=1:n(t)
            tmp_conflict = 0;
            for mm=1:6
                sources = candidate_k_last{t}{j,mm};
                if Ellipse{t-1}{sources(1)}.ind_region ~= Ellipse{t-1}{sources(2)}.ind_region 
                    tmp_conflict = tmp_conflict + fmj{t}(j,mm);
                end
            end
            % 如果邻域内都是单独前景，tmp_conflict可能为0，因此要判断
            if ~isa(tmp_conflict, 'double')
                Fop2 = [ Fop2, tmp_conflict <= 0];
            end
        end
    end
end

%% 可选约束3
if any(use_op_cons==3)
    % %################## 可选择的约束三 #################
    % 融合之后不允许立刻分裂
    for t = s_frame+1:e_frame-1
        for j=1:n(t)
            Fop3 = [ Fop3, sum(fmj{t}(j,:)) + sum(fid{t}(j,:)) <= 1 ];
        end
    end
end

%% 可选约束4是错误的！（已弃用）
    %################## 可选择的约束四 #################（错误的约束，注释掉）
    % 因为单假说前景不是必须被解释，而多假说前景的约束写在矛盾集合处，因此此处注释掉
    %     for t = s_frame+1:e_frame % 2-t 帧要求每个前景内至少有一个被解释（入口解释）
    %         for j=1:n(t)
    %             sum_fid = 0;
    %             sum_fiv = 0;
    %             sum_fij = 0;
    %             % 单独前景必须被解释，即入口和为1
    %             if Ellipse{t}{j}.num_hypoth == 1
    % 
    %                 % sum_fij 为所有迁移到 j 的fij之和（入口和），出口和可以用 sum(fij{t}(j,:)) 表示 
    %                 for ind=1:size(conflict_fij{t-1}{j}, 1)
    %                     sum_fij = sum_fij + fij{t-1}( conflict_fij{t-1}{j}(ind,1), conflict_fij{t-1}{j}(ind,2) );
    %                 end
    % 
    %                 % sum_fid 为所有分裂到包含 j 的 pair 的 fid 之和
    %                 for ind=1:numel(conflict_pair_next_xy{t}{j})/2
    %                     sum_fid = sum_fid + fid{t-1}( conflict_pair_next_xy{t}{j}(ind,1), conflict_pair_next_xy{t}{j}(ind,2) );
    %                     sum_fiv = sum_fiv + fiv{t-1}( conflict_pair_next_xy{t}{j}(ind,1), conflict_pair_next_xy{t}{j}(ind,2) );
    %                 end
    %                 Fop4 = [ Fop4, sum_fij + fsj{t}(j) + sum(fmj{t}(j,:)) + sum_fid + sum_fiv >= 1 ];
    %             else
    %                 % 多目标前景中至少有一个要被解释
    %                 %
    %                 % 这部分写在j到j+jplus为一个前景处 
    %             end
    %         end
    %     end
    % 
    %     t = s_frame; % 同时需要添加第一帧的出口必须被解释
    %     for j=1:n(t)
    %         if Ellipse{t}{j}.num_hypoth == 1
    %             sum_fmj = 0;
    %             % sum_fmj 为所有包含 j 的融合 pair 的 fmj 之和
    %             for ind=1:numel(conflict_pair_last_xy{t}{j})/2
    %                 sum_fmj = sum_fmj + fmj{t+1}( conflict_pair_last_xy{t}{j}(ind,1), conflict_pair_last_xy{t}{j}(ind,2) );
    %             end
    %             Fop4 = [ Fop4, sum(fij{t}(j,:)) + fit{t}(j) + sum(fid{t}(j,:)) + sum(fiv{t}(j,:)) + sum_fmj >= 1 ];
    %         else
    %             % 多目标前景中至少有一个要被解释
    %             %
    %             % 这部分写在j到j+jplus为一个前景处
    %         end
    %     end

%% 可选择约束5 
if any(use_op_cons==5)
    %##################  #################
    % 要求分离split出去的2个子细胞必须在同一前景内
    for t = s_frame:e_frame-1
        for j=1:n(t)
            tmp_conflict = 0;
            for mm=1:6
                sons = candidate_k_next{t}{j,mm};
                if Ellipse{t+1}{sons(1)}.ind_region ~= Ellipse{t+1}{sons(2)}.ind_region % 分裂出去的必须是单独前景
                    tmp_conflict = tmp_conflict + fiv{t}(j,mm);
                end
            end
            % 如果邻域内都是单独前景，tmp_conflict可能为double 0，因此要判断
            if ~isa(tmp_conflict, 'double')
                Fop5 = [ Fop5, tmp_conflict <= 0];
            end
        end
    end
    %#################################################
end
toc
disp('   构建约束3：矛盾约束……');
tic

%% 约束条件3： 矛盾假说集合排除 

conflict = cell(e_frame, 1);
for t = s_frame:e_frame-1
    conflict{t}={};
    j = 1;
    while j<=n(t)    % j遍历
        % 单前景假说则j自增1，继续寻找下个多假说前景
        if Ellipse{t}{j}.num_hypoth==1
            j = j + 1;
            continue;
        end
        
        n_basic = numel(Ellipse{t}{j}.flag_combine);
        if n_basic==1
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
            
        %% A. 此处添加一个多目标前景至少要有一个入口的条件
        if 0
            %##########################################
            % 第一帧的多目标前景出口至少为1（非必要约束）
            if t==s_frame
                sum_fmj = 0;
                sum_fij = 0;
                sum_fit = 0;
                sum_fid = 0;
                sum_fiv = 0;
                for j_muti=j:j+jplus
                    % 求出前景内所有细胞merge出口的和
                    for ind=1:numel(conflict_pair_last_xy{t}{j_muti})/2
                        sum_fmj = sum_fmj + fmj{t+1}( conflict_pair_last_xy{t}{j_muti}(ind,1), conflict_pair_last_xy{t}{j_muti}(ind,2) );
                    end
                    % 求出前景内所有细胞其他 4 种出口的和
                    sum_fij = sum_fij + sum(fij{t}(j_muti,:));
                    sum_fit = sum_fit + fit{t}(j_muti);
                    sum_fid = sum_fid + sum(fid{t}(j_muti,:));
                    sum_fiv = sum_fiv + sum(fiv{t}(j_muti,:));
                end
                Fop6 = [ Fop6, sum_fmj + sum_fij + sum_fit + sum_fid + sum_fiv >= 1 ];

            else
            % 2-t帧的多目标前景入口至少为1（非必要约束）
                sum_fmj = 0;
                sum_fij = 0;
                sum_fsj = 0;
                sum_fid = 0;
                sum_fiv = 0;
                for j_muti=j:j+jplus
                    % 求出前景内所有细胞 split 和 divide 入口的和
                    for ind=1:numel(conflict_pair_next_xy{t}{j_muti})/2
                        sum_fid = sum_fid + fid{t-1}( conflict_pair_next_xy{t}{j_muti}(ind,1), conflict_pair_next_xy{t}{j_muti}(ind,2) );
                        sum_fiv = sum_fiv + fiv{t-1}( conflict_pair_next_xy{t}{j_muti}(ind,1), conflict_pair_next_xy{t}{j_muti}(ind,2) );
                    end
                    % 求出前景内所有细胞其他 3 种入口的和

                    sum_fij_m = 0;
                    % sum_fij 为所有迁移到 j 的fij之和（入口和），出口和可以用 sum(fij{t}(j,:)) 表示 
                    for ind=1:size(conflict_fij{t-1}{j_muti}, 1)
                        sum_fij_m = sum_fij_m + fij{t-1}( conflict_fij{t-1}{j_muti}(ind,1), conflict_fij{t-1}{j_muti}(ind,2) );
                    end

                    % 用 sum_fij_m 入口和代替注释掉的 sum(fij{t-1}(:,j_muti))
                    sum_fij = sum_fij + sum_fij_m;
    %                     sum_fij = sum_fij + sum(fij{t-1}(:,j_muti));
                    sum_fmj = sum_fmj + sum(fmj{t}(j_muti,:));
                    sum_fsj = sum_fsj + fsj{t}(j_muti);
                end
                Fop7 = [ Fop7, sum_fmj + sum_fij + sum_fsj + sum_fid + sum_fiv >= 1 ];
            end
        end
            %##########################################   
            %
            % 以上一段是在一个多目标前景内设置入口至少为1的约束
            % 即保证每个多目标前景都必须被解释，不允许虚景的存在
            %
            %##############################################################

        %% B. conflict表示所有可能的矛盾集合
        use_22_cons = 0;
        if use_22_cons % 新方法，22矛盾约束
            % ----------------------------------------------------------- %
            % 新方法，将矛盾约束写成一个函数，找出前景内的22矛盾，逐一建立条件
            conflict_New = [];
            cc = 1; % 用于计数
            for uu=j:j+jplus-1
                for vv=uu+1:j+jplus
                    if ~isequal( Ellipse{t}{uu}.flag_combine & Ellipse{t}{vv}.flag_combine, zeros(1,n_basic)) % 相与不为0则将其加入矛盾集
                        conflict_New(cc,:) = [uu vv]; % 每一行是一个两两矛盾
                        cc = cc + 1;
                    end
                end
            end
            % 调用外部函数计算这个conflict中的矛盾约束，conflict的行数等于该约束的个数（可以以此验证约束的正确性）
            [ F_j_jlus ] = Calculate_Conflict_Constraint( conflict_New, t, s_frame, fij, fit, fid, fiv, fmj, fsj,...
                conflict_pair_last_xy, conflict_fij, conflict_pair_next_xy);

            F3 = [ F3, F_j_jlus ];

        end
        % ----------------------------------------------------------- %
        if ~use_22_cons % 旧方法有些错误，并不能包含所有的矛盾集合（第一行的for就有问题） 2015.7.7已修复（通过在 CX_Ellipse_Optimal 中加入排序）

            danduan_flag = eye(n_basic); % 单位阵

            for uu=1:n_basic % 遍历单位阵的每一行
                conflict{t}{j,uu} = [];
                for vv=j:j+jplus % 遍历一个前景
                    % 相与不为0则将其加入矛盾集
                    if ~isequal( danduan_flag(uu,:) & Ellipse{t}{vv}.flag_combine, zeros(1,n_basic)) 
                        conflict{t}{j,uu} = [conflict{t}{j,uu}, vv];
                    end
                end

                if numel(conflict{t}{j,uu})==1 % 如果茅盾集中只有1个椭圆，则跳过
                    continue;
                end 
                %##########################此处计算出入口约束###############
                conflict_set = conflict{t}{j,uu}; % 矛盾集合的临时变量
                if t==s_frame % 第一帧考虑出口
                    % ----------------------------------------------- %
                    conflict_sum_out = 0;
                    %### 第一帧一个矛盾集内出口至多为1个
                    for jj=1:numel(conflict_set) % conflict_tmp_set为[21, 23]的形式
                        this_j = conflict_set(jj); % 取到矛盾集中一个细胞的编号
                        sum_fmj = 0;
                        % sum_fmj 为融合出口和
                        tmp_xy = conflict_pair_last_xy{t}{this_j}; % 取到包含当前细胞 ind_j 的pair坐标矩阵
                        for ind=1:numel(tmp_xy)/2
                            sum_fmj = sum_fmj + fmj{t+1}( tmp_xy(ind,1), tmp_xy(ind,2) );
                        end
                        % 单个出口和
                        this_out = sum(fij{t}(this_j,:))+fit{t}(this_j)+sum(fid{t}(this_j,:)) + sum(fiv{t}(this_j,:))+sum_fmj;
                        % 矛盾集全体出口和
                        conflict_sum_out = conflict_sum_out + this_out;
                    end
                    F3 = [ F3, conflict_sum_out <= 1 ];   % 出口唯一
                    % ----------------------------------------------- %
                else
                    % ----------------------------------------------- %
                    % 2—t帧考虑入口
                    conflict_sum_in = 0;
                    %### 一个矛盾集内入口至多为1个
                    for jj=1:numel(conflict_set)
                        this_j = conflict_set(jj); % 取到矛盾集中第一个细胞的编号
                        
                        sum_fid = 0;
                        sum_fiv = 0;
                        sum_fij = 0;
                        % =========================================== %
                        % sum_fij 为所有迁移到 j 的fij之和（入口和），出口和可以用 sum(fij{t}(j,:)) 表示 
                        tmp_xy = conflict_fij{t-1}{this_j};
                        for ind=1:size(tmp_xy, 1)
                            sum_fij = sum_fij + fij{t-1}(tmp_xy(ind,1), tmp_xy(ind,2));
                        end
                        % =========================================== %
                        % sum_fid 为所有分裂出口和
                        tmp_xy = conflict_pair_next_xy{t}{this_j}; % 取到包含当前细胞 ind_j 的pair坐标矩阵
                        for ind=1:size(tmp_xy,1)
                            sum_fid = sum_fid + fid{t-1}( tmp_xy(ind,1), tmp_xy(ind,2) );
                            sum_fiv = sum_fiv + fiv{t-1}( tmp_xy(ind,1), tmp_xy(ind,2) );
                        end
                        % 单个入口和
                        this_in = fsj{t}(this_j)+sum(fmj{t}(this_j,:))+sum_fij+sum_fid+sum_fiv;
                        % 矛盾集全体入口和
                        conflict_sum_in = conflict_sum_in + this_in;
                    end    
                    F4 = [ F4, conflict_sum_in <= 1 ]; %%2—n帧的入口和唯一
                    % ----------------------------------------------- %
                end
                %##########################################################
            end

        end % end if 0
        j = j + jplus + 1; % j直接走出这个前景
            
    end  %## end while
end

toc

%% 导出最终的目标函数与总约束

%#################################
% 各个约束条件的解释如下：
%######### 必要约束 4 个 ##########

% F1：2—t-1 帧中 每个椭圆的 input = output
% F2：2—t   帧中 每个椭圆的 input<=1 && 第 1 帧的 output <=1
% F3：1      帧中 每个矛盾集内 output<=1
% F4：2—t   帧中 每个矛盾集内 inputput<=1

%######### 可选约束 6 个 ##########

% Fop1：2—t   帧中 分裂得到的子细胞不能在同一个前景中 
% Fop2：1—t-1 帧中 不是同一个前景的细胞不能发生融合 
% Fop3：2—t-1 帧中 融合得来的大细胞不允许立刻分裂
% Fop4：1—t   帧中 单假说前景必须被解释（包括 2—t 帧的入口解释和第 1 帧的出口解释）（这个约束是错误的！）
% Fop5：2—t-1 帧中 分离得到的2个细胞必须在同一个前景内 
% Fop6：1      帧中 多假说前景 出口至少为1
% Fop7：2—t   帧中 多假说前景 入口至少为1

%##################################
%##########################################################################
% sdpsettings:      求解器选项设置如下
% http://users.isy.liu.se/johanl/yalmip/pmwiki.php?n=Commands.Sdpsettings

% verbose           打印信息量
% solver            指定求解器
% saveduals         设置为0，不保存对偶变量以节省内存
% usex0             配合 assign 命令给求解器设置初始值，以加快求解速度
% assign            http://users.isy.liu.se/johanl/yalmip/pmwiki.php?n=Commands.Assign
% savesdpafile      保存约束条件和目标函数 http://users.isy.liu.se/johanl/yalmip/pmwiki.php?n=Commands.Savesdpafile
% loadsedumidata    载入约束条件和目标函数 http://users.isy.liu.se/johanl/yalmip/pmwiki.php?n=Commands.Loadsedumidata
      
%##########################################################################

Fbase = [ F1, F2, F3, F4 ];
Foptional = [ Fop1 Fop2 Fop3 Fop5 ]; % 4错误，6、7无需启用，1、2与真实情况有矛盾，因此只考虑3、5（2015.10.3）
% clear Fop1 Fop2 Fop3 Fop4 Fop5 Fop6 Fop7;
Ffull = [ Fbase, Foptional ];


% clear F1 F2 F3 F4 Foptional;
% 保存目标函数与约束条件，清空缓存
% clearvars -except frame Ellipse F object_function fij fid fiv fit fsj fmj s_frame e_frame fai_x_z sum_cost; 

end



function [ F_j_jlus ] = Calculate_Conflict_Constraint( conflict, t, s_frame, fij, fit, fid, fiv, fmj, fsj,...
    conflict_pair_last_xy, conflict_fij, conflict_pair_next_xy)
% 
% 

F_j_jlus = []; % 用来表示j到jplus这一个前景内的矛盾

%##########################此处计算出入口约束###############
if t==s_frame % 第一帧考虑出口
    
    %### 第一帧2个矛盾细胞的出口和至多为1
    for i=1:size(conflict,1)
        sum_fmj = 0;
        j = conflict(i,1); % 取到矛盾集中第一个细胞的编号
        k = conflict(i,2);
        constrain_sum1 = sum( fij{t}(k,:) ) + sum( fij{t}(j,:) ); % 迁移出口和
        constrain_sum2 = fit{t}(k) + fit{t}(j);   % 消失出口和
        constrain_sum3 = sum( fid{t}(k,:) ) + sum( fiv{t}(k,:) ) + sum( fid{t}(j,:) ) + sum( fiv{t}(j,:) );% 分离/分裂出口和
        % sum_fmj 为融合出口和
        tmp_xy_fmj = unique([ conflict_pair_last_xy{t}{j};conflict_pair_last_xy{t}{k} ], 'rows'); % 找出出口矛盾坐标，并去掉重复行 
        for ind=1:size(tmp_xy_fmj,1)
            sum_fmj = sum_fmj + fmj{t+1}( tmp_xy_fmj(ind,1), tmp_xy_fmj(ind,2) );
        end
        F_j_jlus = [ F_j_jlus, constrain_sum1 + constrain_sum2 + constrain_sum3 + sum_fmj <= 1 ];   %%出口唯一
    end
    % 算完return，不需要计算后面
    return;
end

% 2—t帧考虑入口
% 2个矛盾细胞的入口和至多为1
for i=1:size(conflict,1)

    j = conflict(i,1); % 取到矛盾集中第一个细胞的编号
    k = conflict(i,2);
    % ----------------------------------------------- %
    sum_fij = 0;
    % sum_fij 为所有迁移到 j 的fij之和（入口和），出口和可以用 sum(fij{t}(j,:)) 表示 
    tmp_xy_fij =  unique([conflict_fij{t-1}{j};conflict_fij{t-1}{k}], 'rows');
    for ind=1:size(tmp_xy_fij, 1)
        sum_fij = sum_fij + fij{t-1}( tmp_xy_fij(ind,1), tmp_xy_fij(ind,2) );
    end
    % ----------------------------------------------- %
    % 用 sum_fij_m 代替 sum( fij{t-1}(:,ind_j) )
    constrain_sum4 = sum_fij + fsj{t}(j) + sum(fmj{t}(j,:)) + fsj{t}(k) + sum(fmj{t}(k,:)); % 迁移、出现、融合
    % ----------------------------------------------- %
    % sum_fid 为所有分裂出口和                
    sum_fid = 0;
    sum_fiv = 0;
    tmp_xy_fidv = unique([conflict_pair_next_xy{t}{j};conflict_pair_next_xy{t}{k}], 'rows'); % 取到包含当前细胞 ind_j 的pair坐标矩阵
    for ind=1:size(tmp_xy_fidv,1)
        sum_fid = sum_fid + fid{t-1}( tmp_xy_fidv(ind,1), tmp_xy_fidv(ind,2) );
        sum_fiv = sum_fiv + fiv{t-1}( tmp_xy_fidv(ind,1), tmp_xy_fidv(ind,2) );
    end
    % ----------------------------------------------- %
    F_j_jlus = [ F_j_jlus, constrain_sum4 + sum_fid + sum_fiv <= 1 ]; %%2—n帧的入口和唯一
end    

end






