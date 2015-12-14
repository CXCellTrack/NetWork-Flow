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
if 0
    dataset = 'competition';
else
    dataset = 'training';
end
[ ~, trackpath ] = getpath( dataset );

% 使用全局变量
global Fij Fit Fid Fiv Fmj Fsj;
global conflict_fij conflict_pair_last_xy conflict_pair_next_xy n;

load([ trackpath, '\GT\GT_after_hand_tune\GT_Flow_Variables_New.mat']);
load([ trackpath, '\Pair\Pre_data_New.mat'], 'conflict_fij','conflict_pair_last_xy','conflict_pair_next_xy','n','Ellipse','SuperPixel');
frame = numel(Fsj);
frame = 65;

%% 开始检查：主要是入口出口守恒

for t=2:frame-1
    disp(['  检查第',num2str(t),'帧...']);
    for j=1:n(t)
        % 调用函数检查j的入口与出口
        [ eventIn eventOut ] = CX_CheckInOut( t, j );
        
        sum_enter = sum(eventIn);
        sum_leave = sum(eventOut);
        
        % 对每个椭圆的入口出口情况进行判断
        if sum_enter==sum_leave && sum_enter<=1
            % 入口和=出口和，且都为0或都为1，正常跳过
            continue;
        else
            % 其他特殊情况如，如双入口单出口等
            disp(eventIn);
            disp(eventOut);
            
            if exist('SuperPixel','var')
                error(['第', num2str(t), '帧编号为', num2str(j), '的超像素入口为',num2str(sum_enter),'，但出口为',num2str(sum_leave),...
                ';  坐标：',num2str(SuperPixel{t}{j}.centroid),';  BSP: ',num2str(SuperPixel{t}{j}.label)]);
            else
                error(['第', num2str(t), '帧编号为', num2str(j), '的椭圆入口为',num2str(sum_enter),'，但出口为',num2str(sum_leave),...
                    '；  坐标：',num2str([Ellipse{t}{j}.x0, Ellipse{t}{j}.y0])]);
            end
        end
        
    end
end
disp('  需要检查的流程变量满足进出守恒约束');










