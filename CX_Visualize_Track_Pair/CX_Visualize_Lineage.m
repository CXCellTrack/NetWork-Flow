%######################################
%
% 2015.5.30 CX on desk
% 作用：这个函数用于将跟踪完成的细胞上色，便于观察
% 数据存储：上完色的图片保存在‘可视化跟踪标记’
% 依赖关系：调用 CX_RePlot_Ellipse 绘制更新后的椭圆
%          调用 CX_plot 和 CX_plot_event 进行跟踪标记
%###################################### 

%% 绘制可视化跟踪
clear;close all;

% 指定在哪个数据集上进行计算（train or test）
if 1
    dataset = 'competition';
else
    dataset = 'training';
end
[ ~, trackpath ] = getpath( dataset );
trackpath = 'E:\datasets\first_edition\training_datasets\N2DL-HeLa\02_4-16_track';
load([trackpath, '\Pair\Pre_data_New.mat']);

disp('  载入真实流程变量数据...');
load([trackpath, '\GT\GT_after_hand_tune\GT_Flow_Variables_New.mat']);
output_addr = [trackpath, '\GT\'];
s_frame = 1;
e_frame = numel(Fmj);
frame = e_frame - s_frame + 1;

color = colormap(hsv);close('1'); % 选择颜色
nc = numel(color)/3;
color = color(randperm(nc),:);


%% 开始进行循环绘制
im_addr = 'E:\datasets\first_edition\training_datasets\N2DL-HeLa\02\t00.tif'; % 绘图底板
im = imread(im_addr);im = imadjust(im);
im = zeros(size(im));
imshow(im);%imcontrast(gca);
hold on;

for t = s_frame:e_frame
    disp(['  正在处理第 ', num2str(t), ' 帧...']);
    %##### 第一帧需要分配随机彩色 #####
    if t==s_frame
        for j=1:n(t)
            %#####这一段代码反复用到，直接打包#####
            sum_Fmj = 0;
            % sum_fmj 为所有包含 j 的融合 pair 的 fmj 之和
            for ind=1:numel(conflict_pair_last_xy{t}{j})/2
                sum_Fmj = sum_Fmj + Fmj{t+1}( conflict_pair_last_xy{t}{j}(ind,1), conflict_pair_last_xy{t}{j}(ind,2) );
            end
            %####################################
            if sum(Fij{t}(j,:)) + Fit{t}(j) + sum(Fid{t}(j,:)) + sum(Fiv{t}(j,:)) + sum_Fmj ~=0   %%有出口的细胞就分配颜色
                Ellipse{t}{j}.color_index = randi(64);
                % 画个园
                CX_plot_ball( Ellipse{t}{j}, color );
            end
        end
    end
        %###############################
    for j=1:n(t)
        %################# 根据不同的事件把颜色传递给后代 ##################
        %########## 入口事件不出现在第一帧 ##########
        if t~=s_frame
            %########## 出现 ##########
            if Fsj{t}(j)==1   
                indc = randi(nc);
                Ellipse{t}{j}.color_index = indc; 
                CX_plot_ball( Ellipse{t}{j}, color ); % 出现比较特殊，新出现的需要在这里上色
            end
            %########## 融合 ##########
            if sum(Fmj{t}(j,:))==1
                mm = find(Fmj{t}(j,:)==1);
                source = candidate_k_last{t}{j,mm};
                
                % 这是一种异常情况！ 
                Ellipse{t-1}{source(1)} = add_color_index( Ellipse{t-1}{source(1)} );
                Ellipse{t-1}{source(2)} = add_color_index( Ellipse{t-1}{source(2)} );
                
                Ellipse{t}{j}.color_index = [ Ellipse{t-1}{source(1)}.color_index, Ellipse{t-1}{source(2)}.color_index ];
                % 绘制迁移线
                CX_plot_line( Ellipse{t-1}{source(1)}, Ellipse{t}{j}, color );
                CX_plot_line( Ellipse{t-1}{source(2)}, Ellipse{t}{j}, color );
                CX_plot_ball( Ellipse{t}{j}, [1 1 1] );  % 融合得到的需要新画球
            end
        end
        
        %########## 出口事件不出现在最后一帧 ##########
        if t~=e_frame
            %########### 分裂 ##########
            if sum(Fid{t}(j,:))==1    
                mm=find(Fid{t}(j,:)==1);
                % 找到2个子细胞，让它们继承颜色
                son = candidate_k_next{t}{j,mm};
                
                Ellipse{t}{j} = add_color_index(Ellipse{t}{j});
                % 这是一种异常情况！
                if numel(Ellipse{t}{j}.color_index)==2 
                    Ellipse{t+1}{son(1)}.color_index = Ellipse{t}{j}.color_index(1);
                    Ellipse{t+1}{son(2)}.color_index = Ellipse{t}{j}.color_index(2);
                elseif numel(Ellipse{t}{j}.color_index)==1
                    Ellipse{t+1}{son(1)}.color_index = Ellipse{t}{j}.color_index;
                    Ellipse{t+1}{son(2)}.color_index = Ellipse{t}{j}.color_index;
                end

                % 绘制分裂线和新球
                CX_plot_line( Ellipse{t}{j}, Ellipse{t+1}{son(1)}, color );
                CX_plot_ball( Ellipse{t+1}{son(1)}, color );
                CX_plot_line( Ellipse{t}{j}, Ellipse{t+1}{son(2)}, color );
                CX_plot_ball( Ellipse{t+1}{son(2)}, color );
            end

            %########### 分离 ##########
            if sum(Fiv{t}(j,:))==1     
                mm=find(Fiv{t}(j,:)==1);
                % 找到2个分离细胞，让它们继承颜色
                split = candidate_k_next{t}{j,mm};
                
                % 这是一种异常情况！
                Ellipse{t}{j} = add_color_index(Ellipse{t}{j});
                
                if numel(Ellipse{t}{j}.color_index)==2 % 如果找到2种颜色就各自领取，否则就继承
                    Ellipse{t+1}{split(1)}.color_index = Ellipse{t}{j}.color_index(1);
                    Ellipse{t+1}{split(2)}.color_index = Ellipse{t}{j}.color_index(2);
                elseif numel(Ellipse{t}{j}.color_index)==1
                     Ellipse{t+1}{split(1)}.color_index = Ellipse{t}{j}.color_index;
                     Ellipse{t+1}{split(2)}.color_index = Ellipse{t}{j}.color_index;
                end
                % 绘制split线和新球
                CX_plot_line( Ellipse{t}{j}, Ellipse{t+1}{split(1)}, color );
                CX_plot_ball( Ellipse{t+1}{split(1)}, color );
                CX_plot_line( Ellipse{t}{j}, Ellipse{t+1}{split(2)}, color );
                CX_plot_ball( Ellipse{t+1}{split(2)}, color );  
            end
            %########### 消失 ##########
            if Fit{t}(j)==1     
%                 CX_plot_event( Ellipse{t}{j}, 't' );
            end

            %########### 迁移 ##########
            if sum(Fij{t}(j,:))==1
                mm = find(Fij{t}(j,:)==1);
                e_next = candidate_fij{t}(j,mm);
                
                % 这是一种异常情况！
                Ellipse{t}{j} = add_color_index(Ellipse{t}{j});
                Ellipse{t+1}{e_next}.color_index = Ellipse{t}{j}.color_index;
                CX_plot_line( Ellipse{t}{j}, Ellipse{t+1}{e_next}, color );
            end

        end
        
    end

end
hold off;
