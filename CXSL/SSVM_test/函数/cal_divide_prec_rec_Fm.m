function [PRE, REC, FM] = cal_divide_prec_rec_Fm(man_track_txt, res_track_txt)


%% 读入tif图片和track。txt
load(man_track_txt);
load(res_track_txt);
man_dir_path = man_track_txt(1:end-13);
man_dir = dir([man_dir_path,'*.tif']);
res_dir_path = res_track_txt(1:end-13);
res_dir = dir([res_dir_path,'*.tif']);

% 进行数据与处理
man_track = man_track(man_track(:,4)~=0,:);
[~,I] = sort(man_track(:,4));
man_track = man_track(I,:); % 按母细胞编号排序

st = tabulate(man_track(:,4)); % 若发现单个的，则删去
sing = find(st(:,2)==1);
if ~isempty(sing)
    for as=sing'
        man_track(man_track(:,4)==as,:) = [];
    end
end
    
res_track = res_track(res_track(:,4)~=0,:);
[~,I] = sort(res_track(:,4));
res_track = res_track(I,:);

st = tabulate(res_track(:,4)); % 若发现单个的，则删去
if size(st,1)==0
    disp('res_track中没有发现分裂事件！');
    PRE = nan;
    REC = nan;
    FM = nan;
    return
end
sing = find(st(:,2)==1);
if ~isempty(sing)
    for as=sing'
        res_track(res_track(:,4)==as,:) = [];
    end
end
res_track = [res_track, zeros(size(res_track,1),1)];

% truth 数目 predict数目
Tcount = size(man_track,1)/2;
Pcount = size(res_track,1)/2;

%% 统计TP数目
TP = 0;
for h=1:2:size(man_track,1)
    % 求出父子信息
    father = man_track(h,4);
    ftime = man_track(h,2);
    son1 = man_track(h,1);
    son2 = man_track(h+1,1);
    sontime = ftime+1;
    % 在res中找出对应项
    h2list = find(res_track(:,2)==ftime);
    if isempty(h2list)
        continue % 若对应时间的没有，则直接跳过
    end 
    for ii=1:2:numel(h2list)
        h2 = h2list(ii);
        if res_track(h2,5)==1 % 检查之前是否被连接过
            continue
        end
        father2 = res_track(h2, 4);
        ftime2 = res_track(h2, 2);
        son12 = res_track(h2, 1);
        son22 = res_track(h2+1, 1);  
        sontime2 = ftime2+1;
        % 读入2张父亲图片
        pic1 = imread([man_dir_path,man_dir(ftime).name]);
        pic2 = imread([res_dir_path,res_dir(ftime).name]);
        % 接下去判断图的交集
        % --------------------------------------------------------------- %
        loc1 = pic1==father;
        loc2 = pic2==father2;
        if sum(sum(loc1.*loc2))<=0 % 如果2个father交集为0则无效
            continue
        end
        % 读入2张子细胞图片
        pic1 = imread([man_dir_path,man_dir(sontime).name]);
        pic2 = imread([res_dir_path,res_dir(sontime2).name]);
        loc1 = pic1==son1;
        loc2 = pic2==son12;
        if sum(sum(loc1.*loc2))<=0 % 如果2个son1交集为0
            continue
        end
        loc1 = pic1==son2;
        loc2 = pic2==son22;
        if sum(sum(loc1.*loc2))<=0 % 如果2个son2交集为0
            continue
        end
        % --------------------------------------------------------------- %
        % 连接成功，做上标记，以防被后面的继续连接
        fprintf('%d %d %d %d ---> %d %d %d %d\n',man_track(h,:), res_track(h2, 1:4));
        res_track(h2,5) = 1;
        TP = TP + 1;
    end
end

%% 最终结果
PRE = TP/Pcount; fprintf('\n\nprecision:\t%f\n',PRE);
REC = TP/Tcount; fprintf('recall:\t\t%f\n',REC);
FM = 2/(1/PRE+1/REC); fprintf('F-measure:\t%f\n',FM);









