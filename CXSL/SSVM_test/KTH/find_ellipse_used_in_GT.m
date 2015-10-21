function e_used = find_ellipse_used_in_GT()


dataset = 'competition';
[ ~, trackpath ] = getpath( dataset );
% 载入 Ellipse 等数据
load([ trackpath, '\Pair\Pre_data_New.mat']);
% 载入最终的流程变量
track_data_addr = [ trackpath, '\GT\GT_after_hand_tune\GT_Flow_Variables_New.mat'];
load(track_data_addr);

%% 开始检查每个椭圆假说是否参与了流程
e_used  = [];
frame = numel(Fmj);

for t=1:frame
    ind_used = 0;
    for j=1:n(t)

        sum_fij = 0;
        sum_fid = 0;
        sum_fiv = 0;
        sum_fmj = 0;
        sum_enter = 0;
        sum_leave = 0;

        if t>1 % 第一帧没有入口
            % sum_fij 为所有迁移到 j 的fij之和（入口和），出口和可以用 sum(fij{t}(j,:)) 表示 
            for ind=1:size(conflict_fij{t-1}{j}, 1)
                sum_fij = sum_fij + Fij{t-1}( conflict_fij{t-1}{j}(ind,1), conflict_fij{t-1}{j}(ind,2) );
            end

            % sum_fid 为所有分裂到包含 j 的 pair 的 fid 之和（入口）
            for ind=1:numel(conflict_pair_next_xy{t}{j})/2
                sum_fid = sum_fid + Fid{t-1}( conflict_pair_next_xy{t}{j}(ind,1), conflict_pair_next_xy{t}{j}(ind,2) );
                sum_fiv = sum_fiv + Fiv{t-1}( conflict_pair_next_xy{t}{j}(ind,1), conflict_pair_next_xy{t}{j}(ind,2) );
            end
            
            sum_enter = sum_fij + Fsj{t}(j) + sum(Fmj{t}(j,:)) + sum_fid + sum_fiv; % j的总入口和
        end
        
        if t<frame % 最后一帧没有出口
            % sum_fmj 为所有包含 j 的融合 pair 的 fmj 之和（出口）
            for ind=1:numel(conflict_pair_last_xy{t}{j})/2
                sum_fmj = sum_fmj + Fmj{t+1}( conflict_pair_last_xy{t}{j}(ind,1), conflict_pair_last_xy{t}{j}(ind,2) );
            end
            
            sum_leave = sum(Fij{t}(j,:)) + Fit{t}(j) + sum(Fid{t}(j,:)) + sum(Fiv{t}(j,:)) + sum_fmj; % j的总出口和
        end

        
        % 对每个椭圆的入口出口情况进行判断
        if sum_enter || sum_leave
            ind_used = ind_used + 1;
            e_used(t,ind_used) = j; % 这里记录了每帧的椭圆使用情况
        end
        
    end
end










