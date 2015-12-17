%######################################
%
% 2015.5.30 CX on desk
% 作用：这个函数用于将跟踪完成的细胞上色，便于观察
% 数据存储：上完色的图片保存在‘可视化跟踪标记’
% 依赖关系：调用 CX_RePlot_Ellipse 绘制更新后的椭圆
%          调用 CX_fill_color_sp 和 CX_fill_color_event 进行跟踪标记
%###################################### 

%% 绘制可视化跟踪
clear;close all;

% 指定在哪个数据集上进行计算（train or test）
if 1
    dataset = 'competition';
else
    dataset = 'training';
end
[ trainpath, trackpath ] = getpath( dataset );

% 定义全局变量
global Fij Fit Fid Fiv Fmj Fsj;
global conflict_fij conflict_pair_last_xy conflict_pair_next_xy n;

load([trackpath, '\Pair\Pre_data_New.mat']);
 
if 0
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
nc = size(color,1);
rng(1)
color = color(randperm(nc),:);

%% 读入绘图底板
png_addr = [trainpath, '\超像素灰度标签\'];
png_dir = dir([ png_addr, '*.png' ]);
% 设置保存图片的文件夹，先清空其中的fig图片
delete([output_addr,'*.fig']);

%% 开始进行循环绘制
for t = s_frame:e_frame
    pic_name = [ png_addr, png_dir(t).name ];
%     openfig(fig_name, 'new'); %'invisible'
    origin = imread(pic_name);
    im = ones([ size(origin), 3 ]);
    eval_str = '';
    
    disp(['  正在处理', pic_name, '...']);
    %##### 第一帧需要分配随机彩色 #####
    if t==s_frame
        for j=1:n(t)
            % 检查出口
            [ eventIn eventOut ] = CX_CheckInOut( t, j );
            if sum(eventOut)>0
                SuperPixel{t}{j}.color_index = randi(nc); % 有出口的细胞就分配颜色
            end
        end
    end
        
    for j=1:n(t)
        %#################### 上色 #######################
        if isfield( SuperPixel{t}{j}, 'color_index' )
            % 如果椭圆含2个颜色，则说明其为merge后迁移来的，绘作白色
            % 注意，这里不能将merge来的上色，因为运行到这里时还没判断是否为merge
            % merge上色的部分在判断merge处（这里上色的是那些在merge后继续迁移的）
%             assert(numel( SuperPixel{t}{j}.color_index )<=2);
            if numel( SuperPixel{t}{j}.color_index )>=2
                im = CX_fill_color_sp( SuperPixel{t}{j}, im, [1 1 1] ); % 针对merge之后再move的情况
            elseif numel( SuperPixel{t}{j}.color_index )==1
                im = CX_fill_color_sp( SuperPixel{t}{j}, im, color, SuperPixel{t}{j}.color_index );
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
                        eval_str = CX_plot_event_sp( eval_str, SuperPixel{t}{j}, 'dson' );
                    case 4 % 分裂得到的小细胞
                        eval_str = CX_plot_event_sp( eval_str, SuperPixel{t}{j}, 'vson' );
                        
                    case 5 % merge
                        mm = find(Fmj{t}(j,:)==1);
                        source = candidate_k_last{t}{j,mm};
                        SuperPixel{t}{j}.color_index = [ SuperPixel{t-1}{source(1)}.color_index, SuperPixel{t-1}{source(2)}.color_index ];
                        eval_str = CX_plot_event_sp( eval_str, SuperPixel{t}{j}, 'm');   % 融合得来的细胞标记
                        im = CX_fill_color_sp( SuperPixel{t}{j}, im, [1 1 1] );  % 融合得到的需要新上色

                    case 6 % appear
                        eval_str = CX_plot_event_sp( eval_str, SuperPixel{t}{j}, 's');  % 新出现细胞标记
                        SuperPixel{t}{j}.color_index = randi(nc); 
                        im = CX_fill_color_sp( SuperPixel{t}{j}, im, color, SuperPixel{t}{j}.color_index ); % 出现比较特殊，新出现的需要在这里上色
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
                        SuperPixel{t+1}{e_next}.color_index = SuperPixel{t}{j}.color_index;  

                    case 2 % disappear
                        eval_str = CX_plot_event_sp( eval_str, SuperPixel{t}{j}, 't' );

                    case 3 % divide
                        mm=find(Fid{t}(j,:)==1);
                        % 找到2个子细胞，让它们继承颜色（改用标记，颜色采用新的）
                        son = candidate_k_next{t}{j,mm};

                        SuperPixel{t+1}{son(1)}.color_index = randi(nc);
                        SuperPixel{t+1}{son(2)}.color_index = randi(nc);

                        eval_str = CX_plot_event_sp( eval_str, SuperPixel{t}{j}, 'd' );

                    case 4 % split
                        mm=find(Fiv{t}(j,:)==1);
                        % 找到2个分离细胞，让它们继承颜色
                        split = candidate_k_next{t}{j,mm};
                        if numel(SuperPixel{t}{j}.color_index)==2 % 如果找到2种颜色就各自领取，否则（说明无merge就split了，有错误）
                            SuperPixel{t+1}{split(1)}.color_index = SuperPixel{t}{j}.color_index(1);
                            SuperPixel{t+1}{split(2)}.color_index = SuperPixel{t}{j}.color_index(2);
                        else
                            SuperPixel{t+1}{split(1)}.color_index = randi(nc);
                            SuperPixel{t+1}{split(2)}.color_index = randi(nc); % 随机新生
                        end
                        eval_str = CX_plot_event_sp( eval_str, SuperPixel{t}{j}, 'v' );  
                        
                    case 5 % merge
                        eval_str = CX_plot_event_sp( eval_str, SuperPixel{t}{j}, 'mson' );
                end
            end
        end
    end

    % 因为im在不断改变，因此画图语句放入eval_str中，等im不变后再画图（2015.12.14机智！）
    im = drawregionboundaries(origin, im, [0 0 0]); % 画上原始分割线
    imshow(im);hold on;
    eval(eval_str);hold off; 
    
    savename = [output_addr, strrep(png_dir(t).name,'png','fig')];
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

close all


