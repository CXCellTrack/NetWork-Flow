function [ eventIn, eventOut ] = CX_CheckInOut( t, jj )
%
% 2015.12.9 这个文件怎么莫名其妙没下半部分了
% 补上了

% 使用全局变量
global Fij Fit Fid Fiv Fmj Fsj;
global conflict_fij conflict_pair_last_xy conflict_pair_next_xy n;
frame = numel(Fsj);

sum_fij = 0;
sum_fid = 0;
sum_fiv = 0;
sum_fmj = 0;
eventIn = zeros(1,6);
eventOut = zeros(1,6);

if t>1 % 第一帧没有入口
    % sum_fij 为所有迁移到 j 的fij之和（入口和），出口和可以用 sum(fij{t}(j,:)) 表示 
    for ind=1:size(conflict_fij{t-1}{jj}, 1)
        sum_fij = sum_fij + Fij{t-1}( conflict_fij{t-1}{jj}(ind,1), conflict_fij{t-1}{jj}(ind,2) );
    end

    % sum_fid 为所有分裂到包含 j 的 pair 的 fid 之和（入口）
    for ind=1:numel(conflict_pair_next_xy{t}{jj})/2
        sum_fid = sum_fid + Fid{t-1}( conflict_pair_next_xy{t}{jj}(ind,1), conflict_pair_next_xy{t}{jj}(ind,2) );
        sum_fiv = sum_fiv + Fiv{t-1}( conflict_pair_next_xy{t}{jj}(ind,1), conflict_pair_next_xy{t}{jj}(ind,2) );
    end

    eventIn(1) = sum_fij;
    eventIn(2) = 0;
    eventIn(3) = sum_fid;
    eventIn(4) = sum_fiv;
    eventIn(5) = sum(Fmj{t}(jj,:));
    eventIn(6) = Fsj{t}(jj);
    
    sum_enter = sum(eventIn); % j的总入口和
end

if t<frame % 最后一帧没有出口
    % sum_fmj 为所有包含 j 的融合 pair 的 fmj 之和（出口）
    for ind=1:numel(conflict_pair_last_xy{t}{jj})/2
        sum_fmj = sum_fmj + Fmj{t+1}( conflict_pair_last_xy{t}{jj}(ind,1), conflict_pair_last_xy{t}{jj}(ind,2) );
    end

    eventOut(1) = sum(Fij{t}(jj,:));
    eventOut(2) = Fit{t}(jj);
    eventOut(3) = sum(Fid{t}(jj,:));
    eventOut(4) = sum(Fiv{t}(jj,:));
    eventOut(5) = sum_fmj;
    eventOut(6) = 0;
    sum_leave = sum(eventOut); % j的总出口和
end











