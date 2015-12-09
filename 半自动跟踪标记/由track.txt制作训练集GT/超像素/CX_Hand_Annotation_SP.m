%######################################
%
% 2015.12.9 CX on desk
% 作用：这个函数用于从GT TRA中得到灰度标签和我的椭圆编号的对应关系，利用半自动标记 
% 数据存储：将得到的对应矩阵保存为 Label_to_Ellipse.mat 
% 依赖关系：调用 CX_myclick 作为GUI单击的相应函数，可以实现选中/取消
%
%######################################

clear;close all;
[ segpath, trackpath ] = getpath( 'training' );

last = max(strfind(segpath, '\'));
gtpath = [segpath(1:last+2), '_GT\TRA\'];
gt_dir = dir([gtpath, '*.tif']); % gt图片的位置

frame = numel(gt_dir);
frame = 10;

% fig操作标记的图片
fig_path = [trackpath,'\GT\label_and_e\'];
mkdir(fig_path)
fig_dir = dir([fig_path, '\*.fig']);

%% 绘制 TRA 中标签*与椭圆间的位置图
center_gt_path = [trackpath, '\GT\center_gt.mat'];
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
    plot_label2e_SP( stats );
end

%% 载入原始椭圆信息
load([ trackpath, '\Pair\Pre_data_new.mat'], 'SuperPixel','n');

% 载入 label 和 ellipse 的对应关系矩阵
if ~exist([ trackpath, '\GT\Label_to_Ellipse.mat'],'file')
    label2e = cell(frame,1); % label2e 矩阵为label和ellipse的对应关系
else
    load([ trackpath, '\GT\Label_to_Ellipse.mat']);
end

center_sp = cell(frame,1); % 保存SP的中心点
distance = cell(frame,1); % 用来保存SP到*的距离
% 增加标记时，需要对矩阵进行扩大
screen_size = get(0,'ScreenSize');

%% 进行半自动标记
global tmp_label2e global_x global_y last_click;

for t=3:frame % 对frame中有，但目前标记中没有的进行标记
    for j=1:n(t)
        center_sp{t}(j,:) = SuperPixel{t}{j}.centroid;
    end
    %% 单假说前景进行自动标记，每个*找到离自己最近的椭圆
    distance{t} = dist(center_sp{t}, center_gt{t}'); 
    
    for label=1:size(center_gt{t},1) % label为GT标签灰度值
        if isnan( distance{t}(1,label) ) % 如果该列为NaN，则对应为NaN
            label2e{t}(label,1) = NaN;
            continue;   
        end
        % 找出距离*点最近的椭圆编号,因为是距离最近，所以一个*唯一对应一个椭圆；但一个椭圆可能对应2个*
        j = find( distance{t}(:,label) == min(distance{t}(:,label)) ); 
        % 如果*点在椭圆内附近，并且椭圆为单假说前景，则直接对应上
        if distance{t}(j,label)<= stats{t}(label).MajorAxisLength && SuperPixel{t}{j}.num_hypoth ==1
            label2e{t}(label,1) = j; % 第一列存放SP编号
            label2e{t}(label,2) = SuperPixel{t}{j}.label; % 第二列存放SP标签
        else
            label2e{t}(label,1) = 0; % 这个标记为0则说明其需要人工标记
        end
    end
    
    %% 找出未被label对应上的椭圆（包括虚景和多假说前景）进行手动标记
    
    disp(['当前正在处理 ',fig_dir(t).name]);
    openfig( [ trackpath, '\GT\label_and_e\', fig_dir(t).name] );
    % 全屏显示
    set(gcf, 'Position', screen_size);
    
    h = get(gca, 'children');
    % 设置按键事项
    
    tmp_label2e = zeros(50,6); % 这个表用于存放label和e的编号
    global_x = 1;
    global_y = 1;
    last_click = '';
    % 功能为按键->显示
    set(h, 'ButtonDownFcn', @CX_myclick_SP); 
    set(gcf, 'keypressfcn', @pressBlank_SP); % 按空格键负责换行')
    % 需要找出上述椭圆的句柄，才能进行标记
    
    % ------------------------------------------------------------------- %
    % 将没被自动标记上的椭圆绘制为红色，*点绘制为蓝色，以加强视觉区分
    % 1、绘制蓝色点
    label_not_attached = find(label2e{t,1}==0);
    for lab=label_not_attached'
        h_label = findobj(h, 'marker','*','DisplayName',num2str(lab) );
        set(h_label, 'color', 'b');
    end
    % 2、绘制红色圆圈
    bsp = cell2mat( cellfun(@(x) x.label, SuperPixel{t},'un',0)' );
    maxbsp = max(bsp); % 找到这幅图中最大的basic sp（也就是basic sp的总数目）
    sp_not_labed = setdiff( 1:maxbsp, label2e{t}(:,2) );
    for jj=sp_not_labed
        h_e = findobj(h, 'color','g','DisplayName',num2str(jj) );
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
%     label2e{t}(:,2) = []; % 清除label2e的第二列（只在画红色圆时有用）
    for i_row=1:size(tmp_label2e,1)
        row = tmp_label2e(i_row,:);
        tmp = row(row~=0);
        this_sp = tmp(2:end);
        flag = cellfun(@(x) isequal(sort(this_sp), sort(x.label)), SuperPixel{t});
        assert(any(flag==1)) % 断言一定能找到那个label组合
        label2e{t}(row(1),1) = find(flag);
        label2e{t}(row(1),2:numel(this_sp)+1) = this_sp;
    end
    
    if 1
        disp('保存结果到 label2e.mat 中');
        save([ trackpath, '\GT\Label_to_Ellipse.mat'], 'label2e');
    end
    
end
%
% 人工标记使用方法： 红色的椭圆是没被对应上的椭圆，先点击白点，再点击椭圆，就可以将其对应上
%                   如果点错了，可以再点击一下将其取消
%


    
    
    





















        


