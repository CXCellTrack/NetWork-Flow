%######################################
%
% 2015.5.30 CX on desk
% ���ã������������Ԥ������Բ����
% ���ݴ洢�����õ������ݱ���Ϊ Pre_data
% ������ϵ������ CX_ellipse_optimal ������Բ�����Ż�
%
%######################################


% ���� SSVM ѧϰ����֮�󣬾Ͳ���Ҫѧϰ��Щ������ 2015.6.27
clear;close all;

if 1
    dataset = 'competition'; % ѡ��ѵ�����ǲ���
else
    dataset = 'training';
end

%% 1�����������Ϣ

% ����CX-Network ����֮������ĳ�ʼ��Բ����

[ segpath trackpath ] = getpath( dataset );

raw_ellipse_path = [segpath, '\FOI���ͼ2.0\raw_ellipse.mat'];
load( raw_ellipse_path );

frame = numel(ellipse);
% frame = 30;

Ellipse = ellipse(1:frame);
tic;
disp('  ���� CX_Ellipse_Optimal ����...');
% ������蹲��ʱ3�룬ʱ����Ҫ����ȷ���������
Ellipse = CX_Ellipse_Optimal( Ellipse ); % �Ż�֮��õ�����Բ������statusΪ0�ģ����Ѷ���ϳ̶Ȳ�Ľ����˳�
toc;

tic;
% ȥ����cell����ֱ��������Բ����
[ Ellipse, n ] = CX_Cut_Ellipse( Ellipse );
toc;

num_var = zeros(1,frame-1);
for t=1:frame-1
    %           
    % ��Ӧ     divide/split    move       merge     leave   enter
    num_var(t) = n(t)*6*2 +  n(t)*4 + n(t+1)*6 + n(t) + n(t+1);    %%2֮֡��ı�����
end
num_var_sum = sum(num_var);   %%�ܱ�����Ŀ

%% 2��ȷ�� move��divide/split �¼���4���� ����ʱ8.9�룩
tic;
candidate_fij = cell(frame-1,1);

for t=1:frame-1 % ����candidate_fij
%     disp(['����� ',num2str(t),' ֡�� move �¼���4����...']);
    jxy = [];
    kxy = [];
    for j=1:n(t)
        jxy = [ jxy; round([ Ellipse{t}{j}.x0, Ellipse{t}{j}.y0 ]) ];
    end
    for k=1:n(t+1)
        kxy = [ kxy; round([Ellipse{t+1}{k}.x0, Ellipse{t+1}{k}.y0]) ];
    end
        
    distance = dist2( jxy, kxy );
    %############### �����뵹����������,ѡ��Χ�ڵ����4��
    for j=1:n(t)
        [~, ind_k ] = sort(distance(j,:),'ascend');
        % ȡ�����4������Ϊ��ѡ��Ǩ��Ŀ��
        candidate_fij{t}(j,:) = ind_k(1:4);
    end
end

% ---------------------------------------------------------------------- %
% �����Ҫ�޸�predata�����ڴ˴���ͣ�£��޸�candidate_fij����Ϣ %
% ����������������������������������������������������������������������������%
% ---------------------------------------------------------------------- %
candidate_k_next = cell(frame-1,1); %% ��һ֡�еĺ�ѡϸ��
for t=1:frame-1 % ����candidate_k_next
%     disp(['����� ',num2str(t),' ֡�� divide/split �¼���4����...']);
    candidate_k_next{t} = cell(n(t),1);
    for j=1:n(t)
        candidate_k_next{t}{j} = candidate_fij{t}(j,:);
        tmp_combine = combntns(candidate_k_next{t}{j,1}, 2);
        % ��ʱcandidate_k_next{t}{j,1}��Ϊ��4����Բ��ţ����2-7λ����6��pair���
        for tmp_i=2:7
            candidate_k_next{t}{j,tmp_i} = tmp_combine(tmp_i-1,:);
        end
        %######### ����Pij������Բ���������һ����Χ�ڣ��˴���Ϊ5a���������Ϊ0 ##########
    end
    % ֻ����6λ��ϱ�����ɾ����һλ4��Բ���
    candidate_k_next{t}(:,1) = [];
end
toc

%% 3��ȷ�� merge �¼���4���򣨺�ʱ8.9�룩
tic;
candidate_k_last = cell(frame,1); %% ��һ֡�еĺ�ѡϸ��

for t=2:frame
%     disp(['����� ',num2str(t),' ֡�� merge �¼���4����...']);
    candidate_k_last{t} = cell(n(t),1);
    jxy = [];
    ixy = [];
    for j=1:n(t)
        jxy = [ jxy; round([ Ellipse{t}{j}.x0, Ellipse{t}{j}.y0 ]) ];
    end
    for i=1:n(t-1)
        ixy = [ ixy; round([Ellipse{t-1}{i}.x0, Ellipse{t-1}{i}.y0]) ];
    end
        
    distance = dist2( jxy, ixy );
    %############### �����뵹����������,ѡ��Χ�ڵ����4��
    for j=1:n(t)
        [~, ind_i ] = sort(distance(j,:),'ascend');
        candidate_k_last{t}{j,1} = ind_i(1:4);
        tmp_combine = combntns(candidate_k_last{t}{j,1}, 2);
        % 2-7λ����pair���
        for tmp_i=2:7
            candidate_k_last{t}{j,tmp_i} = tmp_combine(tmp_i-1,:);
        end
    end
    % ֻ����6λ��ϱ���
    candidate_k_last{t}(:,1) = [];
end
toc;

%% �����Ϲ��̵ķ�����
% candidate_k_next ����ʵ�ʺ� fid ���ѱ�����Ӧ
% �� fid{t}{i,J}=1 ��ʾ t-1 ʱ�̵�ϸ�� i ���ѵ��� t ʱ�̵�ϸ�� pair J
% �� candidate_k_next{t}{i,J} �п��Ծ��忴�� pair J ��Ӧ��2��ϸ��

% ͬ��candidate_k_last ����ʵ�ʺ� fmj �ںϱ�����Ӧ
% �� fmj{t}{i,J}=1 ��ʾ t ʱ�̵�ϸ�� i ���� t-1 ʱ�̵�ϸ�� pair J �ں϶���
% �� candidate_k_last{t}{i,J} �п��Ծ��忴�� pair J ��Ӧt-1ʱ�̵���2��ϸ��

%% P.S. �� CX_ILP_Pair ��ѡԼ���е� 1��2��5 д�ɸ���Ϊ0��Լ����ʽ��2015.6.15����������ע�� 2015.6.27��

% %################## ��ѡ���Լ��һ&�� #################������ת��Ϊ����Լ����
% % ���ѳ�ȥ�Ĳ�������ͬһ��ǰ����
% for t=1:frame-1
%     for j=1:n(t)
%         for mm=1:6
%             sons = candidate_k_next{t}{j,mm};
%             if Ellipse{t+1}{sons(1)}.ind_region == Ellipse{t+1}{sons(2)}.ind_region % ���ѳ�ȥ�Ĳ�������ͬһ��ǰ����
%                 Pid{t}(j,mm) = 0; % ͨ������Ϊ 0 ����Լ��
%             else
%                 Piv{t}(j,mm) = 0; % �������ͬһ��ǰ���У��򲻿���Ϊsplit����ѡԼ��5��
%             end
%         end
%     end
% end
% 
% %################## ��ѡ���Լ���� #################������ת��Ϊ����Լ����
% % ��������ͬһ��ǰ����ϸ�������ں� ��pair��ʽ��ʵ����������ף�ֻ��Ҫ�Ƚ��ں�pair�е�2ϸ���Ƿ���һ��ǰ���ھͿ�����
% for t=2:frame
%     for j=1:n(t)
%         for mm=1:6
%             sources = candidate_k_last{t}{j,mm};
%             if Ellipse{t-1}{sources(1)}.ind_region ~= Ellipse{t-1}{sources(2)}.ind_region 
%                 Pmj{t}(j,mm) = 0;
%             end
%         end
%     end
% end
% 

%% 4��2-t֡�����ì��pair���� ����/���� ����ʱ11�룩
tic;
% conflict_pair_next{t}{j,1}��ʾt-1ʱ�̷��ѵ�tʱ�̵�pair�а���tʱ��j�����꼯�ϣ������������ì��Լ��
conflict_pair_next_xy = cell(frame,1); 
for t=2:frame
%     disp(['����� ',num2str(t),' ֡�� conflict_pair_next ...']);
    ckn = cell2mat(candidate_k_next{t-1});
    for j=1:n(t)
        [xx, yy] = find(ckn==j);
        yy = ceil(yy/2);
        conflict_pair_next_xy{t}{j,1} = [ xx, yy ];
        % --------------------------------------------------------------- %
    end
end
toc;

%% 5��1��t-1֡�ĳ���ì��pair���� �ϲ�merge ����ʱ11�룩
tic;
% conflict_pair_last{t}{j,1}��ʾtʱ��merge��t+1ʱ�� pair�а���tʱ��j�����꼯�ϣ������������ì��Լ��
conflict_pair_last_xy = cell(frame-1,1); 
for t=1:frame-1
%     disp(['����� ',num2str(t),' ֡�� conflict_pair_last_xy ...']);
    ckl = cell2mat(candidate_k_last{t+1});
    for j=1:n(t)
        [xx, yy] = find(ckl==j);
        yy = ceil(yy/2);
        conflict_pair_last_xy{t}{j,1} = [ xx, yy ];
        % --------------------------------------------------------------- %
    end
end
toc;

%% 6��2��t֡�����ì��fij���ϣ�0.03�룩
tic;
% �� conflict_fij{t} ����ʾ t+1 ʱ�̵�Ǩ�����ì�ܼ���
conflict_fij = cell(frame-1,1); % conflict_fij{t}{j,1}��ʾtʱ���п���Ǩ�Ƶ�t+1ʱ�̵�j������Բ����
for t=1:frame-1
%     disp(['����� ',num2str(t),' ֡�� conflict_fij ...']);
    for next_j=1:n(t+1)   
        % --------------------------------------------------------------- %
        % ע��conflict_fij{t}��һ�ִ�λ�ı�ʾ����Ϊt+1ʱ�̵���Բ��ţ���Ϊtʱ��candidate_fij����Բ����
        % ���ֵ�ʱ��۲����������ɣ��� conflict_fij{1}Ϊ45*1������֪��Ϊt=2��ÿ����Բ��ì�ܼ���
        % conflict_pair_last_xy{1}Ϊ44*1������֪��Ϊt=1��ÿ����Բ�ĳ���ì�ܼ���
        % conflict_pair_next_xy{2}Ϊ45*1������֪��Ϊt=2��ÿ����Բ�����ì�ܼ���
        % --------------------------------------------------------------- %
        [xx, yy] = find( candidate_fij{t}==next_j );
        conflict_fij{t}{next_j,1} = [ xx, yy ];
    end
end
toc;

%% �������ݵ� Pre_data.mat ��
if 1
    disp('  �����ݱ�����PairĿ¼�µ�Pre_data_New��');
    save([trackpath,'\Pair\Pre_data_New.mat'], 'Ellipse', 'candidate_k_last',...
        'candidate_k_next', 'conflict_pair_last_xy', 'conflict_pair_next_xy', 'n', 'num_var', 'num_var_sum',...
        'candidate_fij', 'conflict_fij');
end



