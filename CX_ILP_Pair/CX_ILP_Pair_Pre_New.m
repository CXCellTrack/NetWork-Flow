%######################################
%
% 2015.5.30 CX on desk
% 作用：这个函数用于预处理椭圆数据
% 数据存储：将得到的数据保存为 Pre_data
% 依赖关系：调用 CX_ellipse_optimal 进行椭圆数据优化
%
%######################################


% 采用 SSVM 学习参数之后，就不需要学习这些概率了 2015.6.27
clear;close all;

if 1
    dataset = 'competition'; % 选择训练还是测试
else
    dataset = 'training';
end

%% 1、计算基本信息

% 读入CX-Network 处理之后产生的初始椭圆数据

[ segpath trackpath ] = getpath( dataset );

raw_ellipse_path = [segpath, '\FOI拟合图2.0\raw_ellipse.mat'];
load( raw_ellipse_path );

frame = numel(ellipse);
% frame = 30;

Ellipse = ellipse(1:frame);
tic;
disp('  调用 CX_Ellipse_Optimal 函数...');
% 这个步骤共耗时3秒，时间主要花在确定邻域等上
Ellipse = CX_Ellipse_Optimal( Ellipse ); % 优化之后得到的椭圆不包含status为0的，且已对拟合程度差的进行滤除
toc;

tic;
% 去除空cell，拉直并计算椭圆个数
[ Ellipse, n ] = CX_Cut_Ellipse( Ellipse );
toc;

num_var = zeros(1,frame-1);
for t=1:frame-1
    %           
    % 对应     divide/split    move       merge     leave   enter
    num_var(t) = n(t)*6*2 +  n(t)*4 + n(t+1)*6 + n(t) + n(t+1);    %%2帧之间的变量数
end
num_var_sum = sum(num_var);   %%总变量数目

%% 2、确定 move、divide/split 事件的4邻域 （耗时8.9秒）
tic;
candidate_fij = cell(frame-1,1);

for t=1:frame-1 % 计算candidate_fij
%     disp(['计算第 ',num2str(t),' 帧的 move 事件的4邻域...']);
    jxy = [];
    kxy = [];
    for j=1:n(t)
        jxy = [ jxy; round([ Ellipse{t}{j}.x0, Ellipse{t}{j}.y0 ]) ];
    end
    for k=1:n(t+1)
        kxy = [ kxy; round([Ellipse{t+1}{k}.x0, Ellipse{t+1}{k}.y0]) ];
    end
        
    distance = dist2( jxy, kxy );
    %############### 按距离倒数降序排列,选择范围内的最近4个
    for j=1:n(t)
        [~, ind_k ] = sort(distance(j,:),'ascend');
        % 取最近的4邻域作为候选的迁移目标
        candidate_fij{t}(j,:) = ind_k(1:4);
    end
end

% ---------------------------------------------------------------------- %
% 如果需要修改predata，则在此处暂停下，修改candidate_fij的信息 %
% ！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！%
% ---------------------------------------------------------------------- %
candidate_k_next = cell(frame-1,1); %% 下一帧中的候选细胞
for t=1:frame-1 % 计算candidate_k_next
%     disp(['计算第 ',num2str(t),' 帧的 divide/split 事件的4邻域...']);
    candidate_k_next{t} = cell(n(t),1);
    for j=1:n(t)
        candidate_k_next{t}{j} = candidate_fij{t}(j,:);
        tmp_combine = combntns(candidate_k_next{t}{j,1}, 2);
        % 此时candidate_k_next{t}{j,1}处为那4个椭圆编号，因此2-7位分配6个pair组合
        for tmp_i=2:7
            candidate_k_next{t}{j,tmp_i} = tmp_combine(tmp_i-1,:);
        end
        %######### 计算Pij矩阵，椭圆距离必须在一定范围内（此处设为5a）否则概率为0 ##########
    end
    % 只保留6位组合变量，删除第一位4椭圆编号
    candidate_k_next{t}(:,1) = [];
end
toc

%% 3、确定 merge 事件的4邻域（耗时8.9秒）
tic;
candidate_k_last = cell(frame,1); %% 下一帧中的候选细胞

for t=2:frame
%     disp(['计算第 ',num2str(t),' 帧的 merge 事件的4邻域...']);
    candidate_k_last{t} = cell(n(t),1);
    jxy = [];
    ixy = [];
    for j=1:n(t)
        jxy = [ jxy; round([ Ellipse{t}{j}.x0, Ellipse{t}{j}.y0 ]) ];
    end
    for i=1:n(t-1)
        ixy = [ ixy; round([Ellipse{t-1}{i}.x0, Ellipse{t-1}{i}.y0]) ];
    end
        
    distance = dist2( jxy, ixy );
    %############### 按距离倒数降序排列,选择范围内的最近4个
    for j=1:n(t)
        [~, ind_i ] = sort(distance(j,:),'ascend');
        candidate_k_last{t}{j,1} = ind_i(1:4);
        tmp_combine = combntns(candidate_k_last{t}{j,1}, 2);
        % 2-7位分配pair组合
        for tmp_i=2:7
            candidate_k_last{t}{j,tmp_i} = tmp_combine(tmp_i-1,:);
        end
    end
    % 只保留6位组合变量
    candidate_k_last{t}(:,1) = [];
end
toc;

%% 对以上过程的分析：
% candidate_k_next 矩阵实际和 fid 分裂变量对应
% 即 fid{t}{i,J}=1 表示 t-1 时刻的细胞 i 分裂到了 t 时刻的细胞 pair J
% 在 candidate_k_next{t}{i,J} 中可以具体看到 pair J 对应哪2个细胞

% 同理，candidate_k_last 矩阵实际和 fmj 融合变量对应
% 即 fmj{t}{i,J}=1 表示 t 时刻的细胞 i 是由 t-1 时刻的细胞 pair J 融合而来
% 在 candidate_k_last{t}{i,J} 中可以具体看到 pair J 对应t-1时刻的哪2个细胞

%% P.S. 将 CX_ILP_Pair 可选约束中的 1、2、5 写成概率为0的约束形式（2015.6.15新增）（已注释 2015.6.27）

% %################## 可选择的约束一&五 #################（可以转换为概率约束）
% % 分裂出去的不能仍在同一个前景中
% for t=1:frame-1
%     for j=1:n(t)
%         for mm=1:6
%             sons = candidate_k_next{t}{j,mm};
%             if Ellipse{t+1}{sons(1)}.ind_region == Ellipse{t+1}{sons(2)}.ind_region % 分裂出去的不能仍在同一个前景中
%                 Pid{t}(j,mm) = 0; % 通过概率为 0 进行约束
%             else
%                 Piv{t}(j,mm) = 0; % 如果不在同一个前景中，则不可能为split（可选约束5）
%             end
%         end
%     end
% end
% 
% %################## 可选择的约束二 #################（可以转换为概率约束）
% % 不允许不是同一个前景的细胞发生融合 在pair形式下实现这个很容易，只需要比较融合pair中的2细胞是否在一个前景内就可以了
% for t=2:frame
%     for j=1:n(t)
%         for mm=1:6
%             sources = candidate_k_last{t}{j,mm};
%             if Ellipse{t-1}{sources(1)}.ind_region ~= Ellipse{t-1}{sources(2)}.ind_region 
%                 Pmj{t}(j,mm) = 0;
%             end
%         end
%     end
% end
% 

%% 4、2-t帧的入口矛盾pair集合 分裂/分离 （耗时11秒）
tic;
% conflict_pair_next{t}{j,1}表示t-1时刻分裂到t时刻的pair中包含t时刻j的坐标集合，用来制造入口矛盾约束
conflict_pair_next_xy = cell(frame,1); 
for t=2:frame
%     disp(['计算第 ',num2str(t),' 帧的 conflict_pair_next ...']);
    ckn = cell2mat(candidate_k_next{t-1});
    for j=1:n(t)
        [xx, yy] = find(ckn==j);
        yy = ceil(yy/2);
        conflict_pair_next_xy{t}{j,1} = [ xx, yy ];
        % --------------------------------------------------------------- %
    end
end
toc;

%% 5、1—t-1帧的出口矛盾pair集合 合并merge （耗时11秒）
tic;
% conflict_pair_last{t}{j,1}表示t时刻merge到t+1时刻 pair中包含t时刻j的坐标集合，用来制造出口矛盾约束
conflict_pair_last_xy = cell(frame-1,1); 
for t=1:frame-1
%     disp(['计算第 ',num2str(t),' 帧的 conflict_pair_last_xy ...']);
    ckl = cell2mat(candidate_k_last{t+1});
    for j=1:n(t)
        [xx, yy] = find(ckl==j);
        yy = ceil(yy/2);
        conflict_pair_last_xy{t}{j,1} = [ xx, yy ];
        % --------------------------------------------------------------- %
    end
end
toc;

%% 6、2—t帧的入口矛盾fij集合（0.03秒）
tic;
% 用 conflict_fij{t} 来表示 t+1 时刻的迁移入口矛盾集合
conflict_fij = cell(frame-1,1); % conflict_fij{t}{j,1}表示t时刻有可能迁移到t+1时刻的j处的椭圆坐标
for t=1:frame-1
%     disp(['计算第 ',num2str(t),' 帧的 conflict_fij ...']);
    for next_j=1:n(t+1)   
        % --------------------------------------------------------------- %
        % 注意conflict_fij{t}是一种错位的表示，行为t+1时刻的椭圆编号，列为t时刻candidate_fij的椭圆坐标
        % 区分的时候观察其行数即可，如 conflict_fij{1}为45*1，即可知其为t=2中每个椭圆的矛盾集合
        % conflict_pair_last_xy{1}为44*1，即可知其为t=1中每个椭圆的出口矛盾集合
        % conflict_pair_next_xy{2}为45*1，即可知其为t=2中每个椭圆的入口矛盾集合
        % --------------------------------------------------------------- %
        [xx, yy] = find( candidate_fij{t}==next_j );
        conflict_fij{t}{next_j,1} = [ xx, yy ];
    end
end
toc;

%% 保存数据到 Pre_data.mat 中
if 1
    disp('  将数据保存在Pair目录下的Pre_data_New中');
    save([trackpath,'\Pair\Pre_data_New.mat'], 'Ellipse', 'candidate_k_last',...
        'candidate_k_next', 'conflict_pair_last_xy', 'conflict_pair_next_xy', 'n', 'num_var', 'num_var_sum',...
        'candidate_fij', 'conflict_fij');
end



