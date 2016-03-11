function CXSL_Change_to_AOGMM_without_merge( dataset, flowvars_path )

% 这个不考虑merge及split的复杂情况，merge当作新出现，split当作divide

% 将我得到的跟踪结果转换为cell track challenge上的标准格式
% dataset = 'training';
[ segpath, trackpath ] = getpath( dataset );

% 使用全局变量
global Fij Fit Fid Fiv Fmj Fsj;
global conflict_fij conflict_pair_last_xy conflict_pair_next_xy;

load( flowvars_path );
load([trackpath, '\Pair\Pre_data_New.mat']);
frame = numel(Fmj);

%% 2、将6种事件放在一起，将merge和split替换为move
Final_Table = cell(frame,1); %　Final_Table 是记录各个出口的表，来源记载第三列，分裂而来为1，其余为0

for t=1:frame
    for j=1:numel(Ellipse{t})
        Final_Table{t}(j,1) = 0; % 先置0占位置，后面再进行修改
        [ eventIn, eventOut ] = CX_CheckInOut( t, j );
        % 若j的进出口都为0，说明是废置的假说，用-1进行标记
        if isequal(eventIn,zeros(1,6)) && isequal(eventOut,zeros(1,6))
            Final_Table{t}(j,1) = -1;
        end
        % ================== 检查入口 ================== %
        % ================== 检查出口 ================== %
        ev = find(eventOut);
        if isempty(ev)
            continue
        end
        switch ev
            case 1 % 迁移出去
                nextind = find(Fij{t}(j,:));
                nexte = candidate_fij{t}(j, nextind);
                Final_Table{t}(j,1) = nexte;
            case 2 % ------ 消失 ------ %
                Final_Table{t}(j,1) = 0;
            case 3 % ------ 分裂出去 ------ %
                nextind = find(Fid{t}(j,:));
                sons = candidate_k_next{t}{j, nextind};
                Final_Table{t}(j,1) = sons(1);
                Final_Table{t}(j,2) = sons(2);
                % 要表明其子细胞来源
                Final_Table{t+1}(sons(1),3) = j;
                Final_Table{t+1}(sons(2),3) = j;
            case 4 % ------ 分离出去 ------ %
                nextind = find(Fiv{t}(j,:));
                sons = candidate_k_next{t}{j, nextind};          
                % 若自身是单细胞（说明不是merge而来），则分离无效，修改为分裂
                Final_Table{t}(j,1) = sons(1);
                Final_Table{t}(j,2) = sons(2);
                % 要表明其子细胞来源
                Final_Table{t+1}(sons(1),3) = j;
                Final_Table{t+1}(sons(2),3) = j; 
            case 5 % ------ merge出去 ------ %
                % 这部分放在入口处进行修改
            case 6 % ------ 新出现
                % 对出口无影响
        end
    end
end

%% 3、根据final table产生AOG，导出track_txt
for t=1:frame
    Final_Table{t} = [Final_Table{t}, zeros(size(Final_Table{t},1),3-size(Final_Table{t},2))];
end % 将 Final_Table 统一扩展到3列
       
disp('根据 Final_Table 生成 track_txt...');
track_txt = zeros(500,4);
count = 0;
for t=1:frame
    for j=1:size(Final_Table{t},1)  
        tmpt = t;
        tmpj = j;
        if isfield(Ellipse{t}{j}, 'color') || Final_Table{t}(j,1)== -1
            continue
        end
        fprintf('\n正在跟踪 %d:%d -> ', t,j);
        count = count+1;
        Ellipse{t}{j}.color = count; % 先上色，防止后面的重复
        % -------------------------------------------------------------- %
        lastj = Final_Table{t}(j,3); % 此细胞的祖先
        while 1
            if tmpt==frame % 达到最后一帧，进行记录
                track_txt(count,:) = [count, t-1, tmpt-1, lastj];
                fprintf('=> END ');
                break
            end
            flag = Final_Table{tmpt}(tmpj,1:2)>0;
            if isequal(flag, [1,0])
                tmpj = Final_Table{tmpt}(tmpj,1);
                tmpt = tmpt+1;
                Ellipse{tmpt}{tmpj}.color = count; % 给下一个也上色
                fprintf('%d:%d -> ', tmpt,tmpj);
            else % 除了迁移就是消失和分裂，都属于图走到了叶节点
                track_txt(count,:) = [count, t-1, tmpt-1, lastj];
                if isequal(flag, [0,0])
                    fprintf('Death');
                elseif isequal(flag, [1,1])
                    fprintf('Divide => (%d:%d, %d:%d)',...
                            tmpt+1,Final_Table{tmpt}(tmpj,1), tmpt+1,Final_Table{tmpt}(tmpj,2));
                end
                break
            end
        end
        % -------------------------------------------------------------- %
    end
end
fprintf('\n');
ii = 1;
while ii<=size(track_txt,1)
    if isequal(track_txt(ii,:),[0,0,0,0])
        track_txt(ii,:) = [];
    else
        ii = ii +1;
    end
end

for h=1:size(track_txt,1)
    if track_txt(h,4)
        track_txt(h,4) = Ellipse{track_txt(h,2)}{track_txt(h,4)}.color;
    end
end
       
       
% 保存 track_txt 到文本
filetosave = [trackpath(1:end-11), '_RES\res_track.txt'];
fidin = fopen(filetosave,'wt');
if fidin==-1
    warning('RES目录没有！');
    mkdir([trackpath(1:end-11), '_RES']);
end
for ii=1:size(track_txt,1)
    fprintf(fidin, '%d %d %d %d\n',track_txt(ii,:));
end
fclose(fidin);

%% 3-4、一些辅助测试操作
% ========================================= %
% man_track000 必须为3位数字才行，奇怪
gtpath = [trackpath(1:end-11), '_GT\TRA\'];
gt_dir = dir([gtpath,'*.tif']);
for i=1:numel(gt_dir)
    source = [gtpath,gt_dir(i).name];
    desti = [gtpath,'man_track0',gt_dir(i).name(end-5:end)];
    if strcmp(source,desti)
        break
    end
    movefile(source, desti,'f');
end
gtpath = [trackpath(1:end-11), '_GT\SEG\'];
gt_dir = dir([gtpath,'*.tif']);
for i=1:numel(gt_dir)
    source = [gtpath,gt_dir(i).name];
    desti = [gtpath,'man_seg0',gt_dir(i).name(end-5:end)];
    if strcmp(source,desti)
        break
    end
    movefile(source,desti, 'f');
end
% ========================================= %
% % mask000 也一样
% dirpath = 'E:\datasets\first_edition\training_datasets\N2DL-HeLa\01_RES\';
% gt_dir = dir([dirpath,'*.tif']);
% for i=1:numel(gt_dir)
%     movefile([dirpath,gt_dir(i).name], [dirpath,'mask0',gt_dir(i).name(end-5:end)],'f');
% end 

% % 将gt变为80帧，以便比较
% man_track = 'E:\datasets\first_edition\training_datasets\N2DL-HeLa\01_GT\TRA\man_track.txt';
% file = load(man_track);
% file(file(:,3)>=79, 3) = 79;
% i = 1;
% while i<=size(file,1)
%     if file(i,2)>79
%         file(i,:) = [];
%     else
%         i = i+1;
%     end
% end
% fidin = fopen(man_track,'wt');
% for ii=1:size(file,1)
%     fprintf(fidin, '%d %d %d %d\n',file(ii,:));
% end

%% 4、绘出灰度图像
saveaddr = [trackpath(1:end-11), '_RES\'];
lunkuopath = [segpath, '\FOI提取轮廓\'];
lunkuodir = dir([lunkuopath,'*.tif']);

for t=1:frame
    disp(['绘制图片',num2str(t),'...']);tic
    % ----- 读入原简单分割图 ------- %
    imraw = imread([lunkuopath, lunkuodir(t).name]);
    [height,width] = size(imraw);
    imraw = imfill(imraw,'holes'); % imshow(imraw)
    Label = bwlabeln(imraw);
    % ----------------------------- %
    im = uint16(zeros(height,width));
    
    for j=1:numel(Ellipse{t})  
        if ~isfield(Ellipse{t}{j}, 'color')
            continue
        end
        % 在黑图上画出椭圆
        [im,count] = plot_ellipse_label(im, Label, Ellipse{t}{j}, 1);
        if count==0
            warning([num2str(t),':',num2str(j),' count=0']); 
        end   
    end
    
    name = num2str(t-1);
    if t<=10 % 名字必须为3位
        name = ['0', name];
    end
    imwrite(im, [saveaddr,'mask0',name,'.tif']);toc
end
        
      







            





