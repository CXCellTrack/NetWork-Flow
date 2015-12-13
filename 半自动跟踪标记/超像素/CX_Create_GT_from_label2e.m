%######################################
%
% 2015.12.13 CX on desk
% ���ã�����ű����ڴ� Label_to_Ellipse.mat �еõ����̱�����ʽ��GT
% ���ݴ洢�����õ������̱������󱣴�Ϊ Fid/Fiv/Fij/Fmj/Fsj/Fit
%
%######################################

clear;close all

%% ��ȡ man_track �е���Ϣ

[ segpath, trackpath ] = getpath( 'training' );

last = max(strfind(segpath, '\'));
txtpath = [trackpath(1:last+2), '_GT\TRA\man_track.txt'];
man_track = load( txtpath );
man_track(:,2:3) = man_track(:,2:3) + 1;

% ���� CX_Label_to_Ellipse �������label����Բ�Ķ�Ӧ��ϵ��
load([ trackpath, '\GT\Label_to_Ellipse.mat']);
frame = sum(~isemptycell(label2e));
% ���� pre_data �е���Ϣ
load([ trackpath, '\Pair\Pre_data_New.mat']);
current_tracks = man_track(man_track(:,2)<=frame, :); % starttime��frame�Ժ��Ҫ�ų�

% ��ʼ����������
Fid = cell(frame-1,1);
Fiv = cell(frame-1,1);
Fij = cell(frame-1,1);
Fit = cell(frame-1,1);
for t=1:frame-1
    Fid{t}=zeros(n(t),6);
    Fiv{t}=zeros(n(t),6);
    Fij{t}=zeros(n(t),4);
    Fit{t}=zeros(n(t),1);
end

Fmj = cell(frame,1);
Fsj = cell(frame,1);
for t=2:frame
    Fmj{t}=zeros(n(t),6);
    Fsj{t}=zeros(n(t),1);
end

%% Ѱ��move��Ϣ

move_tracks = current_tracks;
for i_m=1:size(current_tracks,1)
    tmp_row = current_tracks(i_m,:);
    % ����֡���������Χ���ͽ�ȡ��ǰ��֡��
    if tmp_row(3)<frame
        n_m = tmp_row(3);
    else
        n_m = frame;
    end
    % ���м����еĴ���֡��move��Ϊ1
    for m_t=tmp_row(2):n_m-1
        e_last = label2e{m_t}(tmp_row(1),1);
        e_next = label2e{m_t+1}(tmp_row(1),1);
        
        % -------------------------------------------------------- %
        % ��Ϊ����ĳЩ*���޶�Ӧ��Բ�����⣬������Ҫ��������Ϊ��ʧ�¼�
        % ����������Ȼ������ԭʼ�������Ҳֻ������
        if e_last==0 && e_next==0
            % 2�߾�Ϊ0��˵������2֡��û�ж�Ӧ��Բ�������Ա��
        elseif e_last==0 && e_next~=0
            % ǰ��Ϊ0������߳���
            Fsj{m_t+1}(e_next) = 1;
        elseif e_next==0 && e_last~=0
            % ����Ϊ0����ǰ����ʧ
            Fit{m_t}(e_last) = 1;
        elseif e_next~=0 && e_last~=0
            % ---------------------------- %
            % 2�߾���Ϊ0��������Ǩ��
            % ע����Щ���Ǩ�ư�����3�������
            % 1������Ե��㣺move
            % 2������Զ�㣺split
            % 3�����Ե��㣺merge
            % ����2��3��Ǩ����Ҫ�ں���������ɾ�����������ȷ���¼�
            % ---------------------------- %
            % �ҵ� e_next �����Բ�� fij �е��ź� e_next_ind 2015.6.30
            e_next_ind = find( candidate_fij{m_t}(e_last,:)==e_next );
            if isempty(e_next_ind)
                disp(['��',num2str(m_t),'֡��',num2str(e_last),'Ǩ�Ƶ�Ŀ����Բ',num2str(e_next),'����4����ѡ��Բ��']);
                candidate_fij{m_t}(e_last,1) = e_next;
            else
                Fij{m_t}(e_last, e_next_ind) = 1;
            end
        end % 2�߾�Ϊ0����ʲô��������
        % -------------------------------------------------------- %
    end
end

%% Ѱ��appear��Ϣ
app_tracks = current_tracks(current_tracks(:,4)==0, :);
app_tracks = app_tracks(app_tracks(:,2)~=1, :);
for a_t=1:size(app_tracks,1)
    tmp_row = app_tracks(a_t,:);
    % �ҵ��³��ֵ���Բe_app, ���������Ϊ1
    e_app = label2e{tmp_row(2)}(tmp_row(1),1);
    % e_app�п���Ϊ0�����õ��޶�Ӧ��Բʱ���������Ҫ�ж���
    if e_app
        Fsj{tmp_row(2)}(e_app) = 1;
    end
end

%% Ѱ��die��Ϣ
die_tracks = current_tracks(current_tracks(:,3)<frame,:);
tongji = tabulate(current_tracks(:,4));
fathers = tongji(tongji(:,2)==2); % ĸϸ���б�

for h=1:size(die_tracks,1)
    tmp_row = die_tracks(h,:);
    if any(fathers==tmp_row(1)) % ȥ��ĸϸ��
        continue;
    end
    e_die = label2e{tmp_row(3)}(tmp_row(1),1);
    % e_die�п���Ϊ0�����õ��޶�Ӧ��Բʱ���������Ҫ�ж���
    if e_die
        Fit{tmp_row(3)}(e_die) = 1;
    end
end

%% Ѱ��divide��Ϣ

divide_tracks = current_tracks(current_tracks(:,4)~=0, :);
% ����ͳ�Ƴ����ĵ����в�������˳�����У����������
divide_tracks = sortrows(divide_tracks, 4);
tongji = tabulate(divide_tracks(:,4));
not_fa = tongji(tongji(:,2)~=2); % ����һЩ����Ϊ3�����������Ҫɾȥ
ind_delete = [];
for h=1:size(divide_tracks,1)
    if any(not_fa==divide_tracks(h,4))
        ind_delete = [ind_delete, h];
    end
end
divide_tracks(ind_delete,:) = [];


% ���㵱ǰtrack���ܷ��Ѵ���
divide_num = size(divide_tracks,1)/2;
for d_t=1:divide_num
    t = divide_tracks(d_t*2, 2); % tΪ���ѷ�����ʱ��֡+1
    %  �ҳ���ǩ�ϵĸ�����
    label_father = divide_tracks(d_t*2-1, 4);
    label_son1 = divide_tracks(d_t*2-1, 1);
    label_son2 = divide_tracks(d_t*2, 1);
    % ��Ӧ����Բ�ϵĸ�����
    e_father = label2e{t-1}(label_father,1);
    e_son1 = label2e{t}(label_son1,1);
    e_son2 = label2e{t}(label_son2,1);
    % ----------------------------------------- %
    % ����Ƿ�����֮��Ӧ����Բ
    switch num2str(logical([e_father e_son1 e_son2]))
        case '0  1  0' % �޸����ӣ�����Ϊ����
            Fsj{t}(e_son1) = 1;
        case '0  0  1'
            Fsj{t}(e_son2) = 1;
        case '0  1  1'
            Fsj{t}(e_son1) = 1;
            Fsj{t}(e_son2) = 1;
        case '1  0  0' % 2��ϸ����û�Ҽ�������ʧ
            Fit{t-1}(e_father) = 1;
        case '1  1  0' 
            ii = find(candidate_fij{t-1}(e_father,:)==e_son1);
            Fij{t-1}(e_father, ii) = 1;
        case '1  0  1'
            ii = find(candidate_fij{t-1}(e_father,:)==e_son2);
            Fij{t-1}(e_father, ii) = 1;
        case '1  1  1' % Ϊ��������
            % ----------------------------------------- %
            if e_son1==e_son2 % ��Ȼ�����˻�����һ��ǰ����,����Ǩ���¼�
                ii = find(candidate_fij{t-1}(e_father,:)==e_son1);
                Fij{t-1}(e_father, ii) = 1;
                continue;
            end
            real_m =0;
            for mm=1:6
                if isempty( setxor( candidate_k_next{t-1}{e_father, mm}, [e_son1,e_son2] ) )
                    real_m = mm;
                end
            end
            if real_m
                % ���ӷ����¼�
                Fid{t-1}(e_father, real_m) = 1;
            else
                disp(['��',num2str(t-1),'֡����divide��ĸϸ��',num2str(e_father),'��2����ϸ��',num2str([e_son1,e_son2]),'����ĸϸ���ĵ�6����ѡ��Բpair��']);
                candidate_k_next{t-1}{e_father,1} = [e_son1, e_son2];
            end
    end

end

%% ����GT���̱��� 
if 0
    disp('  ����ѵ�����ϵ�GT���̱���...');
    save([ trackpath, '\GT\GT_Flow_Variables_New.mat'], 'Fij','Fid','Fiv','Fit','Fsj','Fmj');
end

% ������ֲ���4�����6�����ڵ���������ʵ��޸������ڵ�����Ҳ���޸�label2e�����ڴ˱���
if 0
    % ע�⣡���޸���candidate����Ҫ��CX_ILP_Pair_Pre_New��������� conflict_pair��conflict_fij
    tic
    disp('���� predata...');
    save([ trackpath, '\Pair\Pre_data_New.mat'], 'Ellipse', 'candidate_k_last',...
        'candidate_k_next', 'conflict_pair_last_xy', 'conflict_pair_next_xy', 'n', 'num_var', 'num_var_sum',...
        'candidate_fij', 'conflict_fij');
    toc
end

    
    
   















