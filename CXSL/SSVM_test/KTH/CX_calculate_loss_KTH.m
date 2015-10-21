% ================================================================== %
%
% CX 2015.10.7
% 这个脚本用于计算KTH得到的精度
% 将KTH中的label与我的椭圆假说对应上，并转化成流程变量的形式
% ================================================================== %


clear;close all

dataset = 'competition';
[ ~, trackpath ] = getpath( dataset );
% 载入predata数据
load([ trackpath, '\Pair\Pre_data_New.mat']); 

last = max(strfind(trackpath, '\'));
KTH_RES_PATH = [ trackpath(1:last+2), '_KTH_RES\'];

% 计算标准答案中使用过的假说（用来标记假说与KTH细胞的对应关系）
tic
e_used = find_ellipse_used_in_GT(); 
toc

% 载入或计算stats（相当于label和e的对应关系）
stats_path = [ KTH_RES_PATH, 'stats.mat' ];
if exist(stats_path, 'file')
    load(stats_path);
else
    % 调用函数进行计算（时间花费较长）
    stats = make_label2e_KTH( KTH_RES_PATH, e_used, Ellipse );
    save(stats_path, 'stats');
end

KTH_track = load([ KTH_RES_PATH, 'res_track.txt']); % 载入KTH的track.txt
KTH_track(:,2:3) = KTH_track(:,2:3) + 1;

%% 分配流程变量
frame = numel(stats);
fij = cell(frame-1,1);
fid = cell(frame-1,1);
fiv = cell(frame-1,1);
fit = cell(frame-1,1);
fsj = cell(frame,1);
fmj = cell(frame,1);

for t = 1:frame-1
    %  t中第m个前景中的第width个细胞   
    fij{t} = zeros(n(t), 4);   %%fij变量矩阵
    fit{t} = zeros(n(t), 1);   %%消失
    fid{t} = zeros(n(t), 6);   %%母细胞
    fiv{t} = zeros(n(t), 6);   %%分裂
end

for t = 2:frame
    fmj{t} = zeros(n(t), 6);   %%融合 ##注意这个矩阵代表（t * t-1），其他都是（t * t+1）
    fsj{t} = zeros(n(t), 1);   %%出现 
end

%% 1、处理move事件
for h=1:size(KTH_track,1)
    row = KTH_track(h,:);
    label = row(1);
    s_frame = row(2);
    e_frame = row(3);
    for t=s_frame:e_frame-1
        e_last = stats{t}(label).e;
        e_next = stats{t+1}(label).e;
        if numel(e_last)==1 && numel(e_next)==1
            ind = find(candidate_fij{t}(e_last,:)==e_next);
            if isempty(ind)
                disp(['  第',num2str(t),'帧 ',num2str(e_last),' 的迁移目标',num2str(e_next),'不在4邻域内！']);
                continue;
            end
            fij{t}(e_last,ind) = 1; % 设置一次迁移事件
        end
    end
end
count_move = sum(KTH_track(:,3) - KTH_track(:,2)); % KTH中迁移出现的次数

%% 2、处理分裂事件
tongji = tabulate(KTH_track(:,4));
tongji = tongji(tongji(:,2)==2); % 这一列是母细胞的label
count_divide = numel(tongji); % KTH中分裂出现的次数
for i_f=1:count_divide
    father = tongji(i_f); % 母细胞label
    sons = find(KTH_track(:,4)==father); % 子细胞label
    t_d = KTH_track(sons(1),2)-1; % 分裂时刻
    
    e_father = stats{t_d}(father).e; % 母细胞假说编号
    e_son1 = stats{t_d+1}(sons(1)).e; % 子细胞假说编号
    e_son2 = stats{t_d+1}(sons(2)).e;
    if isempty(e_father)||isempty(e_son1)||isempty(e_son2)||numel(e_father)>1||numel(e_son1)>1||numel(e_son2)>1
        continue; % 出现各种异常情况，则放弃标记
    end
    
    divide_flag = 0;
    for mm=1:6
        if isempty(setdiff(candidate_k_next{t_d}{e_father,mm}, [e_son1,e_son2]))
            fid{t_d}(e_father,mm) = 1; % 设置一次分裂事件
            divide_flag = 1;
            break;
        end
    end
    if ~divide_flag % 说明6邻域内没找到子细胞pair
        error(['  第',num2str(t_d),'帧 ',num2str(e_father),' 的分裂目标',num2str(e_son1),'和',num2str(e_son2),'不在6邻域内！']);
    end
end

%% 3、处理出现事件和消失事件
app_flag = false(size(KTH_track,1),1);
count_appear = 0;
count_disappear = 0;
for h=1:size(KTH_track,1)
    row = KTH_track(h,:);
        
    if row(2)>1 && ~any(tongji==row(4)) % 不是子细胞且不在第一帧出现，则为出现
        count_appear = count_appear + 1;
        e_app = stats{row(2)}(row(1)).e;
        if numel(e_app)==1
            fsj{row(2)}(e_app) = 1; % 设置出现事件
        end
    end

    if row(3)<frame && ~any(tongji==row(1)) % 如果不是母细胞且提前消失了，则为消失
        count_disappear = count_disappear + 1;
        e_dis = stats{row(3)}(row(1)).e;
        if numel(e_dis)==1
            fit{row(3)}(e_dis) = 1; % 设置出现事件
        end
    end
    
end

%% 载入标准答案进行精度比较
Pcount = zeros(4,1);
Pcount(1) = count_move;
Pcount(2) = count_disappear; % 真实出现的次数在外部计算好
Pcount(3) = count_divide;
Pcount(4) = count_appear;

% 调用函数计算精度（假说精度）
s_frame = 1;
e_frame = numel(stats);
[ PRF COUNT ] = calculate_Loss_Without_Merge( dataset, true, s_frame, e_frame, Pcount, fij, fit, fid, fsj );
% 由于KTH中无merge、split事件，因此调用 calculate_Loss_Without_Merge 来求解精度

if isa(PRF, 'struct')

    disp(PRF.FM) % 只显示F-measure就行了  
    
    % 打印各事件精度
    Preci = struct2array(PRF.Preci)*100;
    Recall = struct2array(PRF.Recall)*100;
    FMeasure = struct2array(PRF.FM)*100;
    
    PRM_for_excel = [ Preci;Recall;FMeasure ];
    COUNT_for_excel = [ COUNT.Tcount'; COUNT.Pcount' ];
    disp(COUNT_for_excel)
end
% ------------------------------------------------------ %
if 0
    txtpath = [ trackpath, '\测试结果记录\KTH.txt'];
    fid = fopen(txtpath, 'w'); fclose(fid);
    save(strrep(txtpath,'txt','mat'), 'PRF','COUNT','fij','fid','fsj','fit'); % 注意修改mat名称
end




