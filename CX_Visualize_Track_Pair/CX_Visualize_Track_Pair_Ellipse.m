%######################################
%
% 2015.5.30 CX on desk
% 作用：这个函数用于将跟踪完成的细胞上色，便于观察
% 数据存储：上完色的图片保存在‘可视化跟踪标记’
% 依赖关系：调用 CX_RePlot_Ellipse 绘制更新后的椭圆
%          调用 CX_fill_color 和 CX_fill_color_event 进行跟踪标记
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

% 定义全局变量
global Fij Fit Fid Fiv Fmj Fsj;
global conflict_fij conflict_pair_last_xy conflict_pair_next_xy n;

load([trackpath, '\Pair\Pre_data_New.mat']);
fig_addr = [trackpath, '\新拟合图\'];
 
if 0
    disp('  载入真实流程变量数据...');
    load([trackpath, '\GT\GT_after_hand_tune\GT_Flow_Variables_New.mat']);
    output_addr = [trackpath, '\GT\GT_after_hand_tune\'];
    s_frame = 1;
    e_frame = numel(Fmj);
else
    disp('  载入经过SSVM learing后得到的分配方案');
    track_data_addr = [trackpath, '\结构化学习\asl_track.mat'];
    output_addr = [trackpath, '\结构化学习\asl_track\'];
    mkdir(output_addr)
    load( track_data_addr );
    % 根据载入的流程变量计算 s_frame 和 e_frame
    s_frame = sum(cellfun(@isempty, Fsj));
    e_frame = numel(Fsj);
end
frame = e_frame - s_frame + 1;

color = colormap(hsv);close('1'); % 选择颜色
nc = size(color,1);
rng(1)
color = color(randperm(nc),:);

%% 读入绘图底板，即新拟合图
fig_dir = dir([ fig_addr, '*.fig' ]);
% 设置保存图片的文件夹，先清空其中的fig图片
delete([output_addr,'*.fig']);
% 即当新拟合图文件夹内的图片数量比目前的frame小时，需要重绘拟合图片
if 0
    CX_RePlot_Ellipse( dataset );
end
% 椭圆是否需要边框，无边框更漂亮，适合作图！有边框用于查错debug
need_contour = true;

%% 开始进行循环绘制
for t = s_frame:e_frame
    fig_name = [ fig_addr, fig_dir(t).name ];
    if need_contour
        openfig(fig_name, 'new'); %'invisible'
    else
        sizeim = [ ]; % 无边框模式下需要填写图片大小
        im = zeros(sizeim);
        imshow(im);
    end
    
    disp(['  正在处理', fig_name, '...']);
    hold on 
    %##### 第一帧需要分配随机彩色 #####
    if t==s_frame
        for j=1:n(t)
            % 检查出口
            [ eventIn eventOut ] = CX_CheckInOut( t, j );
            if sum(eventOut)>0
                Ellipse{t}{j}.color_index = randi(nc); % 有出口的细胞就分配颜色
            end
        end
    end
        
    for j=1:n(t)
        %#################### 上色 #######################
        if isfield( Ellipse{t}{j}, 'color_index' )
            % 如果椭圆含2个颜色，则说明其为merge后迁移来的，绘作白色
            % 注意，这里不能将merge来的上色，因为运行到这里时还没判断是否为merge
            % merge上色的部分在判断merge处（现在已经可以在这里上色了）
            assert(numel( Ellipse{t}{j}.color_index )<=2);
            if numel( Ellipse{t}{j}.color_index )==2
                CX_fill_color( Ellipse{t}{j}, [1 1 1] ); % 针对merge之后再move的情况
            elseif numel( Ellipse{t}{j}.color_index )==1
                CX_fill_color( Ellipse{t}{j}, color, Ellipse{t}{j}.color_index );
            end
        end
    
        %################# 根据不同的事件把颜色传递给后代 ##################
        % 检查入口出口
        [ eventIn eventOut ] = CX_CheckInOut( t, j );
        
        %########## 入口事件不出现在第一帧 ##########
        if t~=s_frame
            ev = find(eventIn);
            if ~isempty(ev) % 加快速度
                switch ev
                    case 3 % 子细胞
                        CX_plot_event( Ellipse{t}{j}, 'dson' );
                    case 4 % 分裂得到的小细胞
                        CX_plot_event( Ellipse{t}{j}, 'vson' );
                        
                    case 5 % merge
                        mm = find(Fmj{t}(j,:)==1);
                        source = candidate_k_last{t}{j,mm};
                        Ellipse{t}{j}.color_index = [ Ellipse{t-1}{source(1)}.color_index, Ellipse{t-1}{source(2)}.color_index ];
                        CX_plot_event( Ellipse{t}{j}, 'm');   % 融合得来的细胞标记
                        CX_fill_color( Ellipse{t}{j}, [1 1 1] );  % 融合得到的需要新上色

                    case 6 % appear
                        CX_plot_event( Ellipse{t}{j}, 's');  % 新出现细胞标记
                        Ellipse{t}{j}.color_index = randi(nc); 
                        CX_fill_color( Ellipse{t}{j}, color, Ellipse{t}{j}.color_index ); % 出现比较特殊，新出现的需要在这里上色
                end
            end
        end
        
        %########## 出口事件不出现在最后一帧 ##########
        if t~=e_frame
            ev = find(eventOut);
            if ~isempty(ev) % 加快速度
                switch ev
                    case 1 % move
                        mm = find(Fij{t}(j,:)==1);
                        e_next = candidate_fij{t}(j,mm);
                        Ellipse{t+1}{e_next}.color_index = Ellipse{t}{j}.color_index;  

                    case 2 % disappear
                        CX_plot_event( Ellipse{t}{j}, 't' );

                    case 3 % divide
                        mm=find(Fid{t}(j,:)==1);
                        % 找到2个子细胞，让它们继承颜色（改用标记，颜色采用新的）
                        son = candidate_k_next{t}{j,mm};

                        Ellipse{t+1}{son(1)}.color_index = randi(nc);
                        Ellipse{t+1}{son(2)}.color_index = randi(nc);

                        CX_plot_event( Ellipse{t}{j}, 'd' );

                    case 4 % split
                        mm=find(Fiv{t}(j,:)==1);
                        % 找到2个分离细胞，让它们继承颜色
                        split = candidate_k_next{t}{j,mm};
                        if numel(Ellipse{t}{j}.color_index)==2 % 如果找到2种颜色就各自领取，否则（说明无merge就split了，有错误）
                            Ellipse{t+1}{split(1)}.color_index = Ellipse{t}{j}.color_index(1);
                            Ellipse{t+1}{split(2)}.color_index = Ellipse{t}{j}.color_index(2);
                        else
                             Ellipse{t+1}{split(1)}.color_index = randi(nc);
                             Ellipse{t+1}{split(2)}.color_index = randi(nc); % 随机新生
                        end
                        CX_plot_event( Ellipse{t}{j}, 'v' );  
                        
                    case 5 % merge
                        CX_plot_event( Ellipse{t}{j}, 'mson' );
                end
            end
        end
    end
    hold off
    savename = [output_addr, strrep(fig_dir(t).name,'tif','fig')];
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




