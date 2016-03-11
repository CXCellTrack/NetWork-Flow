%######################################
%
% 2015.6.6 CX on desk
% 作用：这个函数用于从GT TRA中得到灰度标签和我的椭圆编号的对应关系，利用半自动标记 
% 数据存储：将得到的对应矩阵保存为 Label_to_Ellipse.mat 
% 依赖关系：调用 CX_myclick 作为GUI单击的相应函数，可以实现选中/取消
%
%######################################

clear;close all;
[ ~, trackpath ] = getpath( 'training' );

last = max(strfind(trackpath, '\'));
gtpath = [trackpath(1:last+2), '_GT\TRA\'];

gt_dir = dir([gtpath, '*.tif']);
fig_dir = dir([trackpath,'\新拟合图\*.fig']);
output_dir = dir([trackpath,'\GT\label_and_e\*.fig']);
center_gt_path = [trackpath, '\GT\center_gt.mat'];

frame = numel(fig_dir);

%% 绘制 TRA 中标签*与椭圆间的位置图

if exist(center_gt_path, 'file')
    load(center_gt_path);
else
    center_gt = cell(frame,1);
    stats = cell(frame,1);
    % 先计算label的数据，包括中心坐标和半长轴
    tic
    for t=1:frame
        gt = imread([ gtpath, gt_dir(t).name ]);
        stats{t} = regionprops(gt, 'Centroid', 'MajorAxisLength'); % 计算椭圆的一些信息
        center_gt{t} = cat(1,stats{t}.Centroid);
    end
    toc
    save(center_gt_path, 'stats','center_gt');
end

% ---------------------------------------- %
if 0 % 是否绘制label2e图像
    plot_label2e( stats, fig_dir );
end

%% 载入原始椭圆信息
load([ trackpath, '\Pair\Pre_data_new.mat'], 'Ellipse','n');

% 载入 label 和 ellipse 的对应关系矩阵
if ~exist([ trackpath, '\GT\Label_to_Ellipse.mat'],'file')
    label2e = cell(frame,1); % label2e 矩阵为label和ellipse的对应关系
else
    load([ trackpath, '\GT\Label_to_Ellipse.mat']);
end

center_e = cell(frame,1); % 保存椭圆的中心点
distance = cell(frame,1); % 用来保存椭圆到*的距离
% 增加标记时，需要对矩阵进行扩大
screen_size = get(0,'ScreenSize');

%% 进行半自动标记

for t=92:frame % 对frame中有，但目前标记中没有的进行标记
    for j=1:n(t)
        center_e{t}(j,1) = Ellipse{t}{j}.x0;
        center_e{t}(j,2) = Ellipse{t}{j}.y0;
    end
    
    %% 单假说前景进行自动标记，每个*找到离自己最近的椭圆
    distance{t} = dist(center_e{t}, center_gt{t}'); 
    
    for label=1:size(center_gt{t},1) % label为GT标签灰度值
        if isnan( distance{t}(1,label) ) % 如果该列为NaN，则对应为NaN
            label2e{t}(label,1) = NaN;
            continue;   
        end
        % 找出距离*点最近的椭圆编号,因为是距离最近，所以一个*唯一对应一个椭圆；但一个椭圆可能对应2个*
        j = find( distance{t}(:,label) == min(distance{t}(:,label)) ); 
        % 如果*点在椭圆内附近，并且椭圆为单假说前景，则直接对应上
        if distance{t}(j,label)<= Ellipse{t}{j}.b + stats{t}(label).MajorAxisLength && Ellipse{t}{j}.num_hypoth ==1
            label2e{t}(label,1) = j;
        else
            label2e{t}(label,1) = 0; % 这个标记为0则说明其需要人工标记
        end
    end
    
    %% 找出未被label对应上的椭圆（包括虚景和多假说前景）进行手动标记
    e_not_labed = setdiff( 1:numel(Ellipse{t}), label2e{t}(:,1) );
    label_not_attached = find(label2e{t}==0);
    if isempty(e_not_labed) && isempty(label_not_attached)
        continue;
    end
    
    % 打开当前图片，把上述椭圆标记成红色
    disp(['当前正在处理 ',fig_dir(t).name]);
    openfig( [ trackpath, '\GT\label_and_e\', fig_dir(t).name] );
    set(gcf, 'Position', screen_size);
    % 全屏显示
%     screen_size = get(0,'ScreenSize');
%     set(gcf,'Position',screen_size);
    h = get(gca, 'children');
%     h = findobj(h, 'Type','Line');
    % 设置按键事项
    global tmp_label2e global_x global_y last_click;
    tmp_label2e = zeros(50,2); % 这个表用于存放label和e的编号
    global_x = 1;global_y = 1;last_click = 0;
    % 功能为按键->显示
    set(h, 'ButtonDownFcn', @CX_myclick); 
    % 需要找出上述椭圆的句柄，才能进行标记
    % ------------------------------------------------------------------- %
    % 将没被自动标记上的椭圆绘制为红色，*点绘制为蓝色，以加强视觉区分
    % 1、绘制蓝色点
    for ind = 1:numel(label_not_attached)
        h_label = findobj(h, 'color', 'w', 'DisplayName', num2str(label_not_attached(ind)) );
        set(h_label, 'color', 'b');
    end
    % 2、绘制红色椭圆
    for ind_j=1:numel(e_not_labed)
        h_e = findobj(h, 'color', 'g', 'DisplayName', num2str(e_not_labed(ind_j)) );
        set(h_e, 'color', 'r'); % 红色
    end
    % ------------------------------------------------------------------- %
    
    % 程序暂停，等待处理完图像
    keyin = input('  在控制台上按enter键进入下一帧\n', 's');
    if strcmp(keyin, '')
        close(gcf);
    end

    % 整理标记后的 tmp_label2e 矩阵，去除0和nan
    tmp_label2e = tmp_label2e(tmp_label2e(:,1)~=0, :);
    tmp_label2e = tmp_label2e(~isnan(tmp_label2e(:,1)), :);
    % 将对应关系交给 label2e 矩阵
    label2e{t}(tmp_label2e(:,1),1) = tmp_label2e(:,2);
    clear global;
    
    if 1
        disp('保存结果到 label2e.mat 中');
        save([ trackpath, '\GT\Label_to_Ellipse.mat'], 'label2e');
    end
    
end
%
% 人工标记使用方法： 红色的椭圆是没被对应上的椭圆，先点击白点，再点击椭圆，就可以将其对应上
%                   如果点错了，可以再点击一下将其取消
%


    
    
    





















        


