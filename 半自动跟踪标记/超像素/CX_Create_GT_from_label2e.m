%######################################
%
% 2015.12.13 CX on desk
% 作用：这个脚本用于从 Label_to_Ellipse.mat 中得到流程变量形式的GT
% 数据存储：将得到的流程变量矩阵保存为 Fid/Fiv/Fij/Fmj/Fsj/Fit
%
%######################################

clear;close all

%% 读取 man_track 中的信息

[ segpath, trackpath ] = getpath( 'training' );

last = max(strfind(segpath, '\'));
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
        e_last = label2e{m_t}(tmp_row(1),1);
        e_next = label2e{m_t+1}(tmp_row(1),1);
        
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
    e_app = label2e{tmp_row(2)}(tmp_row(1),1);
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
    e_die = label2e{tmp_row(3)}(tmp_row(1),1);
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
    e_father = label2e{t-1}(label_father,1);
    e_son1 = label2e{t}(label_son1,1);
    e_son2 = label2e{t}(label_son2,1);
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

    
    
   















