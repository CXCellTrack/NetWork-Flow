%######################################
%
% 2015.6.6 CX on desk
% 作用：这个脚本用于从 Label_to_Ellipse.mat 中得到流程变量形式的GT
% 数据存储：将得到的流程变量矩阵保存为 Fid/Fiv/Fij/Fmj/Fsj/Fit
% 依赖关系：调用 
%
%######################################

clear;close all
%% 读取 man_track 中的信息

[ ~, trackpath ] = getpath( 'training' );

last = max(strfind(trackpath, '\'));
txtpath = [trackpath(1:last+2), '_GT\TRA\man_track.txt'];
man_track = load( txtpath );
man_track(:,2:3) = man_track(:,2:3) + 1;

% 读入 CX_Label_to_Ellipse 中求出的label与椭圆的对应关系表
load([ trackpath, '\GT\Label_to_Ellipse.mat']);
frame = sum(~isemptycell(label2e));
% 读入 pre_data 中的信息
load([ trackpath, '\Pair\Pre_data_New.mat']);
current_tracks = man_track(man_track(:,2)<=frame, :); % starttime在frame以后的要排除

% 初始化变量矩阵
Fid = cell(frame-1,1);
Fiv = cell(frame-1,1);
Fij = cell(frame-1,1);
Fit = cell(frame-1,1);
for t=1:frame-1
    Fid{t}=zeros(n(t),6);
    Fiv{t}=zeros(n(t),6);
    Fij{t}=zeros(n(t),4);
    Fit{t}=zeros(n(t),1);
end

Fmj = cell(frame,1);
Fsj = cell(frame,1);
for t=2:frame
    Fmj{t}=zeros(n(t),6);
    Fsj{t}=zeros(n(t),1);
end

%% 寻找move信息

move_tracks = current_tracks;
for i_m=1:size(current_tracks,1)
    tmp_row = current_tracks(i_m,:);
    % 结束帧如果超出范围，就截取当前总帧数
    if tmp_row(3)<frame
        n_m = tmp_row(3);
    else
        n_m = frame;
    end
    % 对中间所有的传递帧，move设为1
    for m_t=tmp_row(2):n_m-1
        e_last = label2e{m_t}(tmp_row(1));
        e_next = label2e{m_t+1}(tmp_row(1));
        
        % -------------------------------------------------------- %
        % 因为存在某些*点无对应椭圆的问题，所以需要将其设置为消失事件
        % 这样设置虽然不满足原始情况，但也只能这样
        if e_last==0 && e_next==0
            % 2者均为0，说明连续2帧都没有对应的圆，不加以标记
        elseif e_last==0 && e_next~=0
            % 前者为0，则后者出现
            Fsj{m_t+1}(e_next) = 1;
        elseif e_next==0 && e_last~=0
            % 后者为0，则前者消失
            Fit{m_t}(e_last) = 1;
        elseif e_next~=0 && e_last~=0
            % ---------------------------- %
            % 2者均不为0，则正常迁移
            % 注意这些点的迁移包含了3种情况：
            % 1、单点对单点：move
            % 2、单点对多点：split
            % 3、多点对单点：merge
            % 其中2和3的迁移需要在后续处理中删除，并添加正确的事件
            % ---------------------------- %
            % 找到 e_next 这个椭圆在 fij 中的排号 e_next_ind 2015.6.30
            e_next_ind = find( candidate_fij{m_t}(e_last,:)==e_next );
            if isempty(e_next_ind)
                disp(['第',num2str(m_t),'帧的',num2str(e_last),'迁移的目标椭圆',num2str(e_next),'不在4个候选椭圆内']);
                candidate_fij{m_t}(e_last,1) = e_next;
            else
                Fij{m_t}(e_last, e_next_ind) = 1;
            end
        end % 2者均为0，则什么都不发生
        % -------------------------------------------------------- %
    end
end

%% 寻找appear信息
app_tracks = current_tracks(current_tracks(:,4)==0, :);
app_tracks = app_tracks(app_tracks(:,2)~=1, :);
for a_t=1:size(app_tracks,1)
    tmp_row = app_tracks(a_t,:);
    % 找到新出现的椭圆e_app, 设置其出现为1
    e_app = label2e{tmp_row(2)}(tmp_row(1));
    % e_app有可能为0（当该点无对应椭圆时），因此需要判断下
    if e_app
        Fsj{tmp_row(2)}(e_app) = 1;
    end
end

%% 寻找die信息
die_tracks = current_tracks(current_tracks(:,3)<frame,:);
tongji = tabulate(current_tracks(:,4));
fathers = tongji(tongji(:,2)==2); % 母细胞列表

for h=1:size(die_tracks,1)
    tmp_row = die_tracks(h,:);
    if any(fathers==tmp_row(1)) % 去除母细胞
        continue;
    end
    e_die = label2e{tmp_row(3)}(tmp_row(1));
    % e_die有可能为0（当该点无对应椭圆时），因此需要判断下
    if e_die
        Fit{tmp_row(3)}(e_die) = 1;
    end
end

%% 寻找divide信息

divide_tracks = current_tracks(current_tracks(:,4)~=0, :);
% 上面统计出来的第四行并不按照顺序排列，因此需排下
divide_tracks = sortrows(divide_tracks, 4);
tongji = tabulate(divide_tracks(:,4));
not_fa = tongji(tongji(:,2)~=2); % 存在一些分裂为3个的情况，需要删去
ind_delete = [];
for h=1:size(divide_tracks,1)
    if any(not_fa==divide_tracks(h,4))
        ind_delete = [ind_delete, h];
    end
end
divide_tracks(ind_delete,:) = [];


% 计算当前track中总分裂次数
divide_num = size(divide_tracks,1)/2;
for d_t=1:divide_num
    t = divide_tracks(d_t*2, 2); % t为分裂发生的时间帧+1
    %  找出标签上的父与子
    label_father = divide_tracks(d_t*2-1, 4);
    label_son1 = divide_tracks(d_t*2-1, 1);
    label_son2 = divide_tracks(d_t*2, 1);
    % 对应到椭圆上的父与子
    e_father = label2e{t-1}(label_father);
    e_son1 = label2e{t}(label_son1);
    e_son2 = label2e{t}(label_son2);
    % ----------------------------------------- %
    % 检查是否有与之对应的椭圆
    switch num2str(logical([e_father e_son1 e_son2]))
        case '0  1  0' % 无父有子，则子为出现
            Fsj{t}(e_son1) = 1;
        case '0  0  1'
            Fsj{t}(e_son2) = 1;
        case '0  1  1'
            Fsj{t}(e_son1) = 1;
            Fsj{t}(e_son2) = 1;
        case '1  0  0' % 2子细胞均没找见，则父消失
            Fit{t-1}(e_father) = 1;
        case '1  1  0' 
            ii = find(candidate_fij{t-1}(e_father,:)==e_son1);
            Fij{t-1}(e_father, ii) = 1;
        case '1  0  1'
            ii = find(candidate_fij{t-1}(e_father,:)==e_son2);
            Fij{t-1}(e_father, ii) = 1;
        case '1  1  1' % 为正常分裂
            % ----------------------------------------- %
            if e_son1==e_son2 % 虽然分裂了还还在一个前景内,视作迁移事件
                ii = find(candidate_fij{t-1}(e_father,:)==e_son1);
                Fij{t-1}(e_father, ii) = 1;
                continue;
            end
            real_m =0;
            for mm=1:6
                if isempty( setxor( candidate_k_next{t-1}{e_father, mm}, [e_son1,e_son2] ) )
                    real_m = mm;
                end
            end
            if real_m
                % 增加分裂事件
                Fid{t-1}(e_father, real_m) = 1;
            else
                disp(['第',num2str(t-1),'帧发生divide的母细胞',num2str(e_father),'的2个子细胞',num2str([e_son1,e_son2]),'不在母细胞的的6个候选椭圆pair内']);
                candidate_k_next{t-1}{e_father,1} = [e_son1, e_son2];
            end
    end

end

%% 寻找merge信息
for t=2:frame-1 
    tongji = tabulate(label2e{t}); % 统计各个椭圆出现的次数
    % 找出发生merge的椭圆集合（可能有多个）
    merge_collect = tongji(tongji(:,2)==2, 1);  % 出现2次的椭圆说明包含了2个*
    merge_collect = setdiff(merge_collect,0);
    % 若无出现2次的，则继续下一轮
    if isempty(merge_collect)
        continue;
    end

    for iind=1:numel( merge_collect ) % 对于集合中的每个椭圆，修改其 Fmj 
        % 找出2个*点
        bigcell = merge_collect(iind);
        label_merge = find(label2e{t} == bigcell)';
        
        % --------------------------------------------------------------- %
        % 判断是否新出现有2种情况
        % 1、新点，即编号很大的点新出现，利用编号是否超过上一帧的最大数判断
        % 2、旧点，编号小，但在上一帧中无对应椭圆，也为新出现（pair2merge中出现nan）
        % 因此需要针对以上2种情况分别判断
        e_son1 = label2e{t-1}(label_merge(1));
        e_son2 = label2e{t-1}(label_merge(2));
        
        if sum([e_son1,e_son2]==0)==1 % 其中一个未标记上，则肯定被判定为出现了
            Fsj{t}(bigcell) = 0; % 需要将出现事件取消掉
            continue;
        end

        flag = isnan([e_son1, e_son2]);
        if isequal(flag, [1 1]) % 2个都是新出现
            Fsj{t}(bigcell) = 0; %?? 暂时还有点问题
        elseif isequal(flag, [0 1]) || isequal(flag, [1 0]) % 有一个为新出现
            Fsj{t}(bigcell) = 0;
        else % 2个源细胞都存在
            % --------------------------------------------------------------- %
            % 如果这2个*点上一帧就在一个椭圆里面，说明这一帧是一个move，无需修改
            if e_son1==e_son2
                continue;
            else
                % 上述2种情况都没有发生，说明的确是merge，则需要查找merge位置
                % 查找merge pair的位置
                real_m = 0;
                for mm=1:6
                    if isempty( setxor([e_son1,e_son2], candidate_k_last{t}{bigcell, mm}) )
                        real_m = mm;
                    end
                end
                if real_m
                    % 修改变量矩阵
                    Fmj{t}( bigcell, real_m) = 1; 
                    % 找到j处细胞在 candidate_fij 中的列号 ind_j
                    j_ind_1 = find( candidate_fij{t-1}(e_son1,:)==bigcell );
                    j_ind_2 = find( candidate_fij{t-1}(e_son2,:)==bigcell );
                    Fij{t-1}(e_son1, j_ind_1) = 0;  % 将之前的迁移行为取消掉
                    Fij{t-1}(e_son2, j_ind_2) = 0;  % 将之前的迁移行为取消掉
                else
                    disp(['第',num2str(t-1),'帧发生merge的源细胞',num2str([e_son1 e_son2]),'不在第',num2str(t),'帧的大细胞',num2str(bigcell),'的6个候选椭圆pair内']);
                    candidate_k_last{t}{bigcell,1} = [e_son1, e_son2];
                end
            end
        end
        
    end
end

%% 寻找split信息
for t=1:frame-1 
    tongji = tabulate(label2e{t}); % 统计各个椭圆出现的次数
    % 找出发生merge的椭圆集合（可能有多个）
    e_merged = tongji(tongji(:,2)==2, 1);  % 出现2次的椭圆说明包含了2个*
    e_merged = setdiff(e_merged, 0);
    num_merged = numel(e_merged); % 找到包含了2个*的椭圆个数（通常只有1个）
    % 如果没有2个重合的则进入下一帧
    if ~num_merged
        continue;
    end
    for i_m=1:num_merged
        
        bigcell = e_merged(i_m); % 临时母椭圆编号
        
        %  2015.7.16 若有2个细胞团直接消失了一个，不满足分裂和分离，则把消失事件删除，视作一次迁移
        if Fit{t}(bigcell)==1 
            Fit{t}(bigcell)=0;
            continue;
        end
        
        label_sources = find(label2e{t}==bigcell)'; % 找到2个小细胞的集合
        e_son1 = label2e{t+1}(label_sources(1));
        e_son2 = label2e{t+1}(label_sources(2));
        % 如果2个label之后对应的不是同一个椭圆，说明split了
        if e_son1 ~= e_son2 
            % 找出分离后的2个椭圆
            real_m =0;
            for mm=1:6
                if isempty( setxor([e_son1, e_son2], candidate_k_next{t}{bigcell, mm}) )
                    real_m = mm;
                end
            end
            if real_m
            % 修改变量矩阵
                Fiv{t}(bigcell, real_m) = 1; 
                % 将原本的迁移事项取消掉
                % 找到子细胞在 candidate_fij 中的列号 ind_j
                son1_ind = find( candidate_fij{t}(bigcell,:)==e_son1 );
                son2_ind = find( candidate_fij{t}(bigcell,:)==e_son2 );
                Fij{t}( bigcell, son1_ind ) = 0;
                Fij{t}( bigcell, son2_ind ) = 0;
            else
                disp(['第',num2str(t),'帧发生split的母细胞',num2str(bigcell),'的2个子细胞',num2str(e_son1),' ',num2str(e_son2),'不在母细胞的的6个候选椭圆pair内']);
                candidate_k_next{t}{bigcell,1} = [e_son1, e_son2];
            end
            
        end
        % 
        % 如果2个label之后对应的还是同一个椭圆，说明move了
        % 这个部分已经在merge处处理过了，此处无需操作
        %
    end
end

%% 保存GT流程变量 
if 0
    disp('  保存训练集合的GT流程变量...');
    save([ trackpath, '\GT\GT_Flow_Variables_New.mat'], 'Fij','Fid','Fiv','Fit','Fsj','Fmj');
end

% 如果出现不在4邻域或6邻域内的情况，可适当修改邻域内的数，也可修改label2e，并在此保存
if 0
    % 注意！：修改了candidate后，需要到CX_ILP_Pair_Pre_New中重新求解 conflict_pair和conflict_fij
    tic
    disp('保存 predata...');
    save([ trackpath, '\Pair\Pre_data_New.mat'], 'Ellipse', 'candidate_k_last',...
        'candidate_k_next', 'conflict_pair_last_xy', 'conflict_pair_next_xy', 'n', 'num_var', 'num_var_sum',...
        'candidate_fij', 'conflict_fij');
    toc
end

    
    
   















