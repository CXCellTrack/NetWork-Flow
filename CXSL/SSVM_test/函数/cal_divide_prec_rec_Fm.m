function [PRE, REC, FM] = cal_divide_prec_rec_Fm(man_track_txt, res_track_txt)


%% ����tifͼƬ��track��txt
load(man_track_txt);
load(res_track_txt);
man_dir_path = man_track_txt(1:end-13);
man_dir = dir([man_dir_path,'*.tif']);
res_dir_path = res_track_txt(1:end-13);
res_dir = dir([res_dir_path,'*.tif']);

% ���������봦��
man_track = man_track(man_track(:,4)~=0,:);
[~,I] = sort(man_track(:,4));
man_track = man_track(I,:); % ��ĸϸ���������

st = tabulate(man_track(:,4)); % �����ֵ����ģ���ɾȥ
sing = find(st(:,2)==1);
if ~isempty(sing)
    for as=sing'
        man_track(man_track(:,4)==as,:) = [];
    end
end
    
res_track = res_track(res_track(:,4)~=0,:);
[~,I] = sort(res_track(:,4));
res_track = res_track(I,:);

st = tabulate(res_track(:,4)); % �����ֵ����ģ���ɾȥ
if size(st,1)==0
    disp('res_track��û�з��ַ����¼���');
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

% truth ��Ŀ predict��Ŀ
Tcount = size(man_track,1)/2;
Pcount = size(res_track,1)/2;

%% ͳ��TP��Ŀ
TP = 0;
for h=1:2:size(man_track,1)
    % ���������Ϣ
    father = man_track(h,4);
    ftime = man_track(h,2);
    son1 = man_track(h,1);
    son2 = man_track(h+1,1);
    sontime = ftime+1;
    % ��res���ҳ���Ӧ��
    h2list = find(res_track(:,2)==ftime);
    if isempty(h2list)
        continue % ����Ӧʱ���û�У���ֱ������
    end 
    for ii=1:2:numel(h2list)
        h2 = h2list(ii);
        if res_track(h2,5)==1 % ���֮ǰ�Ƿ����ӹ�
            continue
        end
        father2 = res_track(h2, 4);
        ftime2 = res_track(h2, 2);
        son12 = res_track(h2, 1);
        son22 = res_track(h2+1, 1);  
        sontime2 = ftime2+1;
        % ����2�Ÿ���ͼƬ
        pic1 = imread([man_dir_path,man_dir(ftime).name]);
        pic2 = imread([res_dir_path,res_dir(ftime).name]);
        % ����ȥ�ж�ͼ�Ľ���
        % --------------------------------------------------------------- %
        loc1 = pic1==father;
        loc2 = pic2==father2;
        if sum(sum(loc1.*loc2))<=0 % ���2��father����Ϊ0����Ч
            continue
        end
        % ����2����ϸ��ͼƬ
        pic1 = imread([man_dir_path,man_dir(sontime).name]);
        pic2 = imread([res_dir_path,res_dir(sontime2).name]);
        loc1 = pic1==son1;
        loc2 = pic2==son12;
        if sum(sum(loc1.*loc2))<=0 % ���2��son1����Ϊ0
            continue
        end
        loc1 = pic1==son2;
        loc2 = pic2==son22;
        if sum(sum(loc1.*loc2))<=0 % ���2��son2����Ϊ0
            continue
        end
        % --------------------------------------------------------------- %
        % ���ӳɹ������ϱ�ǣ��Է�������ļ�������
        fprintf('%d %d %d %d ---> %d %d %d %d\n',man_track(h,:), res_track(h2, 1:4));
        res_track(h2,5) = 1;
        TP = TP + 1;
    end
end

%% ���ս��
PRE = TP/Pcount; fprintf('\n\nprecision:\t%f\n',PRE);
REC = TP/Tcount; fprintf('recall:\t\t%f\n',REC);
FM = 2/(1/PRE+1/REC); fprintf('F-measure:\t%f\n',FM);









