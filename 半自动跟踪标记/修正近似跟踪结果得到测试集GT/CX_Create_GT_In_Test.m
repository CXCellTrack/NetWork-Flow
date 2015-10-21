% ================================================================== %
%
% CX 2015.7.21
% 这个脚本用于根据 CX_Hand_Annotation_In_Test 得到的结果（包括删除和添加）
% 将一个较为接近 gt 的 trackdata 修正为 ground truth
% 主要步骤包括：1 删去错误轨迹
%              2 增加新的正确轨迹
%
% ================================================================== %

clear;close all;

%% 载入测试集合上的跟踪结果
if 1
    dataset = 'competition';
else
    dataset = 'training';
end

[ ~, trackpath ] = getpath( dataset );

load([ trackpath, '\Pair\Pre_data_New.mat']);
track_data_addr = [ trackpath, '\结构化学习\Tracking_Data.mat'];
if strcmp(dataset,'training') % 修改训练集答案地址在这里
    track_data_addr = [ trackpath, '\GT\GT_Flow_Variables_New.mat'];
end

load( track_data_addr );
frame = numel(Fmj);
% frame = 20;
% 复制载入的结果
Fij_c = Fij;
Fit_c = Fit;
Fid_c = Fid;
Fiv_c = Fiv;
Fmj_c = Fmj;
Fsj_c = Fsj;

%% 打开绘制完成的跟踪图片，制作GT（在 CX_Hand_Annotation_In_Test 中进行操作）
load([ trackpath, '\GT\Hand_GT_New.mat']); 
GT_delete = GT_delete_s;
GT_move = GT_move_s;

%% 1：增加或更改新的轨迹，首先要删除1的去路和2的来路，再添加1->2的轨迹
for t=1:numel(GT_move)
    % 为空则跳过，加快速度
    if isempty(GT_move{t})
        continue;
    end
    
    for ind=1:size(GT_move{t},1)
        % 找到一行记录的内容
        rowdata = GT_move{t}(ind,:);
        if isequal(rowdata, [0 0 0 0]) % 这一行为空记录则跳过
            break;
        end

        %% A 删除要修改的轨迹2端的去路和来路
        % ---------------- 删去j的去路 -------------- %
        for ii=1:2
            j = GT_move{t}(ind,ii); % j为GT_move一行中的前2个
            if j % j=0时为k出现，j无去路
                % 1 迁移出口置0
                Fij_c{t}(j,:) = zeros(1,4);
                % 2 divide/split出口置0
                Fid_c{t}(j,:) = zeros(1,6);
                Fiv_c{t}(j,:) = zeros(1,6);
                % 3 消失出口置0
                Fit_c{t}(j) = 0;
                % 4 merge出口置0
                uu = conflict_pair_last_xy{t}{j};
                for i=1:size(uu,1)
                    Fmj_c{t+1}(uu(i,1),uu(i,2)) = 0;
                end
                disp([ '第', num2str(t), '帧编号为', num2str(j), '的椭圆去路已删除！']);
            end
        end
        % ---------------- 删去k的来路 -------------- %
        for ii=3:4
            k = GT_move{t}(ind,ii); % k为GT_move一行中的后2个
            if k % k=0时为j消失，k无来路
                % 1 迁移入口置0
                uu = conflict_fij{t}{k};
                for i=1:size(uu,1)
                    Fij_c{t}(uu(i,1),uu(i,2)) = 0;
                end
                % 2 merge入口置0
                Fmj_c{t+1}(k,:) = zeros(1,6);
                % 3 divide/split入口置0
                uu = conflict_pair_next_xy{t+1}{k};
                for i=1:size(uu,1)
                    Fid_c{t}(uu(i,1),uu(i,2)) = 0;
                    Fiv_c{t}(uu(i,1),uu(i,2)) = 0;
                end
                % 4 出现入口置0
                Fsj_c{t+1}(k) = 0;
                disp([ '第', num2str(t+1), '帧编号为', num2str(k), '的椭圆来路已删除！']);
            end
        end
        % ------------------------------------------- %
        
        %% B 添加正确的轨迹
        % =============== 添加新的轨迹 =============== %
        % 判断事件
        j1 = rowdata(1);
        j2 = rowdata(2);
        k1 = rowdata(3);
        k2 = rowdata(4);
        
        disp('  增加新的轨迹...');
        
        abcd = num2str(rowdata~=0);
        switch abcd
            case '1  0  0  0'
                Fit_c{t}(j1) = 1;
                disp([ '  第', num2str(t), '帧的', num2str(j1), '发生disappear事件']);
            case '0  0  1  0'
                Fsj_c{t+1}(k1) = 1;
                disp([ '  第', num2str(t+1), '帧的', num2str(k1), '发生appear事件']);
            case '1  0  1  0'
                mm = find(candidate_fij{t}(j1,:)==k1);
                if isempty(mm)
                    msgbox([ '  第', num2str(t), '帧的', num2str(j1), '迁移的目标', num2str(k1), '不在四邻域内！']);
                    candidate_fij{t}(j1,1) = k1; % 修改邻域
                end
                Fij_c{t}(j1,mm) = 1;   
                disp([ '  第', num2str(t), '帧的', num2str(j1), '迁移为', '第', num2str(t+1), '帧的', num2str(k1)]);
            case '1  0  1  1'
                % ------------------------------------------------------- %
                % divide/split需要找出子细胞pair位置,都先记作divide，若有split再手动更改
                mmtrue = 0;
                for mm=1:6
                    if isempty( setdiff(candidate_k_next{t}{j1,mm}, [k1 k2]) )
                        mmtrue = mm;
                    end
                end
                % 判断是否在四邻域内找到了[k1 k2]
                if mmtrue
%                     if Ellipse{t+1}{k1}.ind_region~=Ellipse{t+1}{k2}.ind_region % 不再同一区域为分裂，否则为分离
                        Fid_c{t}(j1,mmtrue) = 1;
                        disp([ '  第', num2str(t), '帧的', num2str(j1), '分裂为',...
                            '第', num2str(t+1), '帧的', num2str(k1), '和', num2str(k2)]);
%                     else
%                         Fiv_c{t}(j1,mmtrue) = 1;
%                         disp([ '  第', num2str(t), '帧的', num2str(j1), '分离为',...
%                             '第', num2str(t+1), '帧的', num2str(k1), '和', num2str(k2)]);
%                     end        
                else
                    msgbox(['  GT_move的第',num2str(t),'帧第',num2str(ind),'行的divide/split标记无法在候选椭圆中找到！'])
                    candidate_k_next{t}{j1,1} = [k1 k2];
                end
                % ------------------------------------------------------- %
            case '1  1  1  0'
                % ------------------------------------------------------- %
                % merge 需要先找出在源细胞pair位置
                mmtrue = 0;
                for mm=1:6
                    if isempty( setdiff(candidate_k_last{t+1}{k1,mm}, [j1 j2]) )
                        mmtrue = mm;
                    end
                end
                % 判断是否在四邻域内找到了[j1 j2]
                if mmtrue
                    Fmj_c{t+1}(k1,mmtrue) = 1;
                    disp([ '  第', num2str(t), '帧的', num2str(j1), '和', num2str(j2), '合并为',...
                        '第', num2str(t+1), '帧的', num2str(k1)]);
                else
                    error(['  GT_move的第',num2str(t),'帧第',num2str(ind),'行的merge标记无法在候选椭圆中找到！'])
                end
                % ------------------------------------------------------- %
            otherwise
                msgbox(['  GT_move的第',num2str(t),'帧第',num2str(ind),'行记录不是标准的细胞事件！']);
                candidate_k_last{t+1}{k1,1} = [j1 j2];
        end
        
        disp(' ');
        % ============================================ %
        
    end
end
     
%% 2：删除标记为删除的椭圆的来路和去路（这些是错误标记，且无法修正为正确的）
for t=1:numel(GT_delete)    
    % 为空则跳过，加快速度
    if isempty(GT_delete{t})
        continue;
    end
    
    for ind=1:numel(GT_delete{t})
        j = GT_delete{t}(ind);
        if ~j % j=0则跳过
            continue;
        end

        % 第一种情况：删除来路
        if t>1
            % 1 迁移入口置0
            uu = conflict_fij{t-1}{j};
            for i=1:size(uu,1)
                Fij_c{t-1}(uu(i,1),uu(i,2)) = 0;
            end
            % 2 merge入口置0
            Fmj_c{t}(j,:) = zeros(1,6);
            % 3 divide/split入口置0
            uu = conflict_pair_next_xy{t}{j};
            for i=1:size(uu,1)
                Fid_c{t-1}(uu(i,1),uu(i,2)) = 0;
                Fiv_c{t-1}(uu(i,1),uu(i,2)) = 0;
            end
            % 4 出现入口置0
            Fsj_c{t}(j) = 0;
        end
        disp([ '第', num2str(t), '帧编号为', num2str(j), '的椭圆来路已删除！']);

        % 第二种情况：删除去路
        if t<frame
            % 1 迁移出口置0
            Fij_c{t}(j,:) = zeros(1,4);
            % 2 divide/split出口置0
            Fid_c{t}(j,:) = zeros(1,6);
            Fiv_c{t}(j,:) = zeros(1,6);
            % 3 消失出口置0
            Fit_c{t}(j) = 0;
            % 4 merge出口置0
            uu = conflict_pair_last_xy{t}{j};
            for i=1:size(uu,1)
                Fmj_c{t+1}(uu(i,1),uu(i,2)) = 0;
            end
        end
        disp([ '第', num2str(t), '帧编号为', num2str(j), '的椭圆去路已删除！']);
    end
end
disp(' ');
disp('========================================');
% -----------------------------

%% 保存制作完成的GT
%  为了不改变后面计算精度的代码，此处把结果保存为Fij的形式
Fij = Fij_c;
Fit = Fit_c;
Fid = Fid_c;
Fiv = Fiv_c;
Fmj = Fmj_c;
Fsj = Fsj_c;

% 变量窗内手动添加divide/split/merge事件，再进行保存
% 目前有。。。
% 2015.7.22已经可以自动配置，无需逐条记录，只需记录split即可
% ---- split ---- %
% t  i  j1  j2  mm
% 64 75 


if 1
    disp('保存人工修改过后的GT流程变量');
    save([ trackpath, '\GT\GT_after_hand_tune\GT_Flow_Variables_New.mat'],...
        'Fij','Fit','Fid','Fiv','Fmj','Fsj');
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

