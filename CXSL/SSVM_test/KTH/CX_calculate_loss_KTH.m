% ================================================================== %
%
% CX 2015.10.7
% ����ű����ڼ���KTH�õ��ľ���
% ��KTH�е�label���ҵ���Բ��˵��Ӧ�ϣ���ת�������̱�������ʽ
% ================================================================== %


clear;close all

dataset = 'competition';
[ ~, trackpath ] = getpath( dataset );
% ����predata����
load([ trackpath, '\Pair\Pre_data_New.mat']); 

last = max(strfind(trackpath, '\'));
KTH_RES_PATH = [ trackpath(1:last+2), '_KTH_RES\'];

% �����׼����ʹ�ù��ļ�˵��������Ǽ�˵��KTHϸ���Ķ�Ӧ��ϵ��
tic
e_used = find_ellipse_used_in_GT(); 
toc

% ��������stats���൱��label��e�Ķ�Ӧ��ϵ��
stats_path = [ KTH_RES_PATH, 'stats.mat' ];
if exist(stats_path, 'file')
    load(stats_path);
else
    % ���ú������м��㣨ʱ�仨�ѽϳ���
    stats = make_label2e_KTH( KTH_RES_PATH, e_used, Ellipse );
    save(stats_path, 'stats');
end

KTH_track = load([ KTH_RES_PATH, 'res_track.txt']); % ����KTH��track.txt
KTH_track(:,2:3) = KTH_track(:,2:3) + 1;

%% �������̱���
frame = numel(stats);
fij = cell(frame-1,1);
fid = cell(frame-1,1);
fiv = cell(frame-1,1);
fit = cell(frame-1,1);
fsj = cell(frame,1);
fmj = cell(frame,1);

for t = 1:frame-1
    %  t�е�m��ǰ���еĵ�width��ϸ��   
    fij{t} = zeros(n(t), 4);   %%fij��������
    fit{t} = zeros(n(t), 1);   %%��ʧ
    fid{t} = zeros(n(t), 6);   %%ĸϸ��
    fiv{t} = zeros(n(t), 6);   %%����
end

for t = 2:frame
    fmj{t} = zeros(n(t), 6);   %%�ں� ##ע������������t * t-1�����������ǣ�t * t+1��
    fsj{t} = zeros(n(t), 1);   %%���� 
end

%% 1������move�¼�
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
                disp(['  ��',num2str(t),'֡ ',num2str(e_last),' ��Ǩ��Ŀ��',num2str(e_next),'����4�����ڣ�']);
                continue;
            end
            fij{t}(e_last,ind) = 1; % ����һ��Ǩ���¼�
        end
    end
end
count_move = sum(KTH_track(:,3) - KTH_track(:,2)); % KTH��Ǩ�Ƴ��ֵĴ���

%% 2����������¼�
tongji = tabulate(KTH_track(:,4));
tongji = tongji(tongji(:,2)==2); % ��һ����ĸϸ����label
count_divide = numel(tongji); % KTH�з��ѳ��ֵĴ���
for i_f=1:count_divide
    father = tongji(i_f); % ĸϸ��label
    sons = find(KTH_track(:,4)==father); % ��ϸ��label
    t_d = KTH_track(sons(1),2)-1; % ����ʱ��
    
    e_father = stats{t_d}(father).e; % ĸϸ����˵���
    e_son1 = stats{t_d+1}(sons(1)).e; % ��ϸ����˵���
    e_son2 = stats{t_d+1}(sons(2)).e;
    if isempty(e_father)||isempty(e_son1)||isempty(e_son2)||numel(e_father)>1||numel(e_son1)>1||numel(e_son2)>1
        continue; % ���ָ����쳣�������������
    end
    
    divide_flag = 0;
    for mm=1:6
        if isempty(setdiff(candidate_k_next{t_d}{e_father,mm}, [e_son1,e_son2]))
            fid{t_d}(e_father,mm) = 1; % ����һ�η����¼�
            divide_flag = 1;
            break;
        end
    end
    if ~divide_flag % ˵��6������û�ҵ���ϸ��pair
        error(['  ��',num2str(t_d),'֡ ',num2str(e_father),' �ķ���Ŀ��',num2str(e_son1),'��',num2str(e_son2),'����6�����ڣ�']);
    end
end

%% 3����������¼�����ʧ�¼�
app_flag = false(size(KTH_track,1),1);
count_appear = 0;
count_disappear = 0;
for h=1:size(KTH_track,1)
    row = KTH_track(h,:);
        
    if row(2)>1 && ~any(tongji==row(4)) % ������ϸ���Ҳ��ڵ�һ֡���֣���Ϊ����
        count_appear = count_appear + 1;
        e_app = stats{row(2)}(row(1)).e;
        if numel(e_app)==1
            fsj{row(2)}(e_app) = 1; % ���ó����¼�
        end
    end

    if row(3)<frame && ~any(tongji==row(1)) % �������ĸϸ������ǰ��ʧ�ˣ���Ϊ��ʧ
        count_disappear = count_disappear + 1;
        e_dis = stats{row(3)}(row(1)).e;
        if numel(e_dis)==1
            fit{row(3)}(e_dis) = 1; % ���ó����¼�
        end
    end
    
end

%% �����׼�𰸽��о��ȱȽ�
Pcount = zeros(4,1);
Pcount(1) = count_move;
Pcount(2) = count_disappear; % ��ʵ���ֵĴ������ⲿ�����
Pcount(3) = count_divide;
Pcount(4) = count_appear;

% ���ú������㾫�ȣ���˵���ȣ�
s_frame = 1;
e_frame = numel(stats);
[ PRF COUNT ] = calculate_Loss_Without_Merge( dataset, true, s_frame, e_frame, Pcount, fij, fit, fid, fsj );
% ����KTH����merge��split�¼�����˵��� calculate_Loss_Without_Merge ����⾫��

if isa(PRF, 'struct')

    disp(PRF.FM) % ֻ��ʾF-measure������  
    
    % ��ӡ���¼�����
    Preci = struct2array(PRF.Preci)*100;
    Recall = struct2array(PRF.Recall)*100;
    FMeasure = struct2array(PRF.FM)*100;
    
    PRM_for_excel = [ Preci;Recall;FMeasure ];
    COUNT_for_excel = [ COUNT.Tcount'; COUNT.Pcount' ];
    disp(COUNT_for_excel)
end
% ------------------------------------------------------ %
if 0
    txtpath = [ trackpath, '\���Խ����¼\KTH.txt'];
    fid = fopen(txtpath, 'w'); fclose(fid);
    save(strrep(txtpath,'txt','mat'), 'PRF','COUNT','fij','fid','fsj','fit'); % ע���޸�mat����
end




