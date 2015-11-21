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

load([trackpath, '\Pair\Pre_data_New.mat']);
fig_addr = [trackpath, '\新拟合图\'];
 
if 1
    disp('  载入真实流程变量数据...');
    load([trackpath, '\GT\GT_after_hand_tune\GT_Flow_Variables_New.mat']);
    output_addr = [trackpath, '\GT\GT_after_hand_tune\'];
    s_frame = 1;
    e_frame = numel(Fmj);
else
    disp('  载入经过SSVM learing后得到的分配方案');
    track_data_addr = [trackpath, '\结构化学习\Tracking_Data.mat'];
    output_addr = [trackpath, '\Pair\可视化跟踪标记\'];
    mkdir(output_addr)
    load( track_data_addr );
    % 根据载入的流程变量计算 s_frame 和 e_frame
    s_frame = sum(cellfun(@isempty, Fsj));
    e_frame = numel(Fsj);
end
frame = e_frame - s_frame + 1;

color = colormap(hsv);close('1'); % 选择颜色
color_numbel = numel(color)/3;
color = color(randperm(color_numbel),:);

% 读入绘图底板，即新拟合图
fig_dir = dir([ fig_addr, '*.fig' ]);
% 设置保存图片的文件夹，先清空其中的fig图片
delete([output_addr,'*.fig']);
% 即当新拟合图文件夹内的图片数量比目前的frame小时，需要重绘拟合图片
if 0
    CX_RePlot_Ellipse( dataset );
end

%% 开始进行循环绘制
for t = s_frame:e_frame
    fig_name = [ fig_addr, fig_dir(t).name ];
    openfig(fig_name, 'new'); %'invisible'
    disp(['  正在处理', fig_name, '...']);
    hold;
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
                end
            end
        end
        %###############################
    for j=1:n(t)
        %#################### 上色 #######################
        if isfield( Ellipse{t}{j}, 'color_index' )
            % 如果椭圆含2个颜色，则说明其为merge后迁移来的，绘作白色
            % 注意，这里不能将merge来的上色，因为运行到这里时还没判断是否为merge
            % merge上色的部分在判断merge处
            if numel( Ellipse{t}{j}.color_index )==2
                CX_plot( Ellipse{t}{j}, [1 1 1] );
            elseif numel( Ellipse{t}{j}.color_index )==1
                CX_plot( Ellipse{t}{j}, color(Ellipse{t}{j}.color_index,:) );
            end
        end
    
        %################# 根据不同的事件把颜色传递给后代 ##################
        
        %########## 入口事件不出现在第一帧 ##########
        if t~=s_frame
            %########## 出现 ##########
            if Fsj{t}(j)==1   
                indc = randi(color_numbel);
                CX_plot_event( Ellipse{t}{j}, 's');  % 新出现细胞标记
                Ellipse{t}{j}.color_index = indc; 
                CX_plot( Ellipse{t}{j}, color(Ellipse{t}{j}.color_index,:) ); % 出现比较特殊，新出现的需要在这里上色
            end
            %########## 融合 ##########
            if sum(Fmj{t}(j,:))==1
                mm = find(Fmj{t}(j,:)==1);
                source = candidate_k_last{t}{j,mm};
                Ellipse{t}{j}.color_index = [ Ellipse{t-1}{source(1)}.color_index, Ellipse{t-1}{source(2)}.color_index ];
                CX_plot_event( Ellipse{t}{j}, 'm');   % 融合得来的细胞标记
                CX_plot( Ellipse{t}{j}, [1 1 1] );  % 融合得到的需要新上色
            end
        end
        
        %########## 出口事件不出现在最后一帧 ##########
        if t~=e_frame
            %########### 分裂 ##########
            if sum(Fid{t}(j,:))==1    
                mm=find(Fid{t}(j,:)==1);
                % 找到2个子细胞，让它们继承颜色
                son = candidate_k_next{t}{j,mm};
                
                Ellipse{t+1}{son(1)}.color_index = Ellipse{t}{j}.color_index;
                Ellipse{t+1}{son(2)}.color_index = Ellipse{t}{j}.color_index;

                CX_plot_event( Ellipse{t}{j}, 'd' );
            end

            %########### 分离 ##########
            if sum(Fiv{t}(j,:))==1     
                mm=find(Fiv{t}(j,:)==1);
                % 找到2个分离细胞，让它们继承颜色
                split = candidate_k_next{t}{j,mm};
                if numel(Ellipse{t}{j}.color_index)==2 % 如果找到2种颜色就各自领取，否则就继承
                    Ellipse{t+1}{split(1)}.color_index = Ellipse{t}{j}.color_index(1);
                    Ellipse{t+1}{split(2)}.color_index = Ellipse{t}{j}.color_index(2);
                else
                     Ellipse{t+1}{split(1)}.color_index = Ellipse{t}{j}.color_index;
                     Ellipse{t+1}{split(2)}.color_index = Ellipse{t}{j}.color_index;
                end
                CX_plot_event( Ellipse{t}{j}, 'v' );  
            end
            %########### 消失 ##########
            if Fit{t}(j)==1     
                CX_plot_event( Ellipse{t}{j}, 't' );
            end

            %########### 迁移 ##########
            if sum(Fij{t}(j,:))==1
                mm = find(Fij{t}(j,:)==1);
                e_next = candidate_fij{t}(j,mm);
                Ellipse{t+1}{e_next}.color_index = Ellipse{t}{j}.color_index;
            end

        end
        
    end
    hold;
    savename = [output_addr, fig_dir(t).name];
    saveas(1,savename);
    close('1');
end

%% 图片显示
if 0
    screen_size = get(0,'ScreenSize');
    colored_fig_dir = dir([ output_addr, '*.fig' ]);  
    for t = 1:20
        fig_name = [ output_addr, colored_fig_dir(t).name ];
        openfig(fig_name, 'new', 'visible');
        set(gcf, 'Position', screen_size);
    end
end