% ======================================================================= %
% 2015.8.27
% 在得到一组流程变量以后：Fij Fid Fit Fiv Fsj Fmj
% 1、如果流程变量是由cplex求解得到的，则正常情况下满足各种约束
% 2、如果流程变量是由人工标记再经过相关处理得到的，则不一定满足相关约束（主要是入口出口守恒约束）
%
% 本脚本的作用就是检查得到的流程变量是否满足必要的约束条件
% 不满足的话则调用 CX_Visualize_Track_Pair_New 进行上色时会报错，因此需要在此处检查
% ======================================================================= %

clear;close all;
% 载入GT流程变量（主要检查这个）  载入pre_data
if 1
    dataset = 'competition';
else
    dataset = 'training';
end
[ ~, trackpath ] = getpath( dataset );

load([ trackpath, '\GT\GT_after_hand_tune\GT_Flow_Variables_New.mat']);
load([ trackpath, '\Pair\Pre_data_New.mat']);
frame = numel(Fmj);
% frame = 40;

%% 开始检查：主要是入口出口守恒

for t=2:frame-1
    disp(['  检查第',num2str(t),'帧...']);
    for j=1:n(t)
        sum_fij = 0;
        sum_fid = 0;
        sum_fiv = 0;
        sum_fmj = 0;

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
            
            eventIn(1) = sum_fij;
            eventIn(2) = Fsj{t}(j);
            eventIn(3) = sum(Fmj{t}(j,:));
            eventIn(4) = sum_fid;
            eventIn(5) = sum_fiv;
            sum_enter = sum_fij + Fsj{t}(j) + sum(Fmj{t}(j,:)) + sum_fid + sum_fiv; % j的总入口和
        end
        
        if t<frame % 最后一帧没有出口
            % sum_fmj 为所有包含 j 的融合 pair 的 fmj 之和（出口）
            for ind=1:numel(conflict_pair_last_xy{t}{j})/2
                sum_fmj = sum_fmj + Fmj{t+1}( conflict_pair_last_xy{t}{j}(ind,1), conflict_pair_last_xy{t}{j}(ind,2) );
            end
            
            eventOut(1) = sum(Fij{t}(j,:));
            eventOut(2) = Fit{t}(j);
            eventOut(3) = sum(Fid{t}(j,:));
            eventOut(4) = sum(Fiv{t}(j,:));
            eventOut(5) = sum_fmj;

            sum_leave = sum(Fij{t}(j,:)) + Fit{t}(j) + sum(Fid{t}(j,:)) + sum(Fiv{t}(j,:)) + sum_fmj; % j的总出口和
        end
        

        % 对每个椭圆的入口出口情况进行判断
        if sum_enter==sum_leave && sum_enter<=1
            % 入口和=出口和，且都为0或都为1，正常跳过
            continue;
        else
            % 其他特殊情况如，如双入口单出口等
            disp(eventIn);
            disp(eventOut);
            error(['第', num2str(t), '帧编号为', num2str(j), '的椭圆入口为',num2str(sum_enter),'，但出口为',num2str(sum_leave),...
                '；  坐标：',num2str([Ellipse{t}{j}.x0, Ellipse{t}{j}.y0])]);
        end
        
    end
end
disp('  需要检查的流程变量满足进出守恒约束');










