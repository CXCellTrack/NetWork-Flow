function [ cost_for_train PRF COUNT ] = CX_Calculate_Loss( dataset, exist_GT, s_frame, e_frame, fij, fit, fid, fiv, fmj, fsj )

% ���㾫��
%% ��������

[ ~, trackpath ] = getpath( dataset );

load([ trackpath, '\Pair\Pre_data_New.mat'], 'n','conflict_pair_last_xy','conflict_pair_next_xy','conflict_fij');
gt_flow_addr = [ trackpath, '\GT\GT_after_hand_tune\GT_Flow_Variables_New.mat'];

% �����׼��GT������Ҫ��GT���ڵ�����£�ʹ��bool exist_GT����ʶGT�Ƿ���ڣ�
if exist_GT
    load( gt_flow_addr );
else
    cost_for_train = [];
    PRF = [];
    COUNT = [];
    return;
end

%% 4.������ʧ����
% ======================================================================= %
% �������¼�����ʧ���������ݹ�ʽ��10�����ⲿ��Ҫ����Ŀ�꺯����
% ======================================================================= %
% ============== delta(f,f*) ============== %
TP.fij = 0;
TP.fit = 0;
TP.fid = 0;
TP.fiv = 0;
TP.fmj = 0;
TP.fsj = 0;

% PN������Խ����TF����GT���
%    T  F
% P  TP FP
% N  FN TN
% 
% precision = tp/(tp+fp)
% recall = tp/(tp+fn);

Tcount = zeros(6,1); % tp+tn
Pcount = zeros(6,1); % tp+fn
% ��ʧ���������� false negative ���㣬�� f*=1 && f=0 ʱ������ʧ
for t = s_frame:e_frame-1
    TP.fij = TP.fij + sum(sum( Fij{t}.*fij{t} ));
    Tcount(1) = Tcount(1) + sum(sum( Fij{t})); % ͳ�Ƴ� Fij{t}��1�ĸ���
    Pcount(1) = Pcount(1) + sum(sum( fij{t})); % �����г��ֵĸ���
    
    TP.fit = TP.fit + sum(sum( Fit{t}.*fit{t} ));
    Tcount(2) = Tcount(2) + sum(sum( Fit{t}));
    Pcount(2) = Pcount(2) + sum(sum( fit{t}));
    
    TP.fid = TP.fid + sum(sum( Fid{t}.*fid{t} ));
    Tcount(3) = Tcount(3) + sum(sum( Fid{t}));
    Pcount(3) = Pcount(3) + sum(sum( fid{t}));
    
    TP.fiv = TP.fiv + sum(sum( Fiv{t}.*fiv{t} ));
    Tcount(4) = Tcount(4) + sum(sum( Fiv{t}));
    Pcount(4) = Pcount(4) + sum(sum( fiv{t}));
end
for t = s_frame+1:e_frame
    TP.fmj = TP.fmj + sum(sum( Fmj{t}.*fmj{t} ));
    Tcount(5) = Tcount(5) + sum(sum( Fmj{t}));
    Pcount(5) = Pcount(5) + sum(sum( fmj{t}));
    
    TP.fsj = TP.fsj + sum(sum( Fsj{t}.*fsj{t} ));
    Tcount(6) = Tcount(6) + sum(sum( Fsj{t}));
    Pcount(6) = Pcount(6) + sum(sum( fsj{t}));
end

% ============= ע�⣺�����еķ����ǰ����е�count�������ڼ��� ============== %
event_TP = TP.fij + TP.fit + TP.fid + TP.fiv + TP.fmj + TP.fsj;
 
%% ���㣨GT�в����õļ�˵ ����ȴ�����ã���ɵ����
% ======================================================================= %
% �����龰����ì�ܼ�˵������������� f* ȫΪ0����f�г�����1����Ҫ�������
% ======================================================================= %

% ����2�����֡����ڣ�����Բj��ʵ���Ϊ0���ҷ������Ϊ1�������һ����ʧ
fd_TP = 0;
fd_Pcount = 0; % ͳ�Ʋ����С��龰�����ִ���
fd_Tcount = 0; % ͳ��GT�С��龰��������ڳ���Ϊ0����Բ�����ֵĴ���

for t = s_frame+1:e_frame
    for j=1:n(t)  
        % ������ڱ���
        sum_fid = 0;
        sum_fiv = 0;
        % ��ʵ��ڱ���
        sum_Fid = 0;
        sum_Fiv = 0;
        % sum_fid Ϊ���з��ѵ����� j �� pair �� fid ֮��
        for ind=1:numel(conflict_pair_next_xy{t}{j})/2
            sum_fid = sum_fid + fid{t-1}( conflict_pair_next_xy{t}{j}(ind,1), conflict_pair_next_xy{t}{j}(ind,2) );
            sum_fiv = sum_fiv + fiv{t-1}( conflict_pair_next_xy{t}{j}(ind,1), conflict_pair_next_xy{t}{j}(ind,2) );
            
            sum_Fid = sum_Fid + Fid{t-1}( conflict_pair_next_xy{t}{j}(ind,1), conflict_pair_next_xy{t}{j}(ind,2) );
            sum_Fiv = sum_Fiv + Fiv{t-1}( conflict_pair_next_xy{t}{j}(ind,1), conflict_pair_next_xy{t}{j}(ind,2) );
        end
        % ����Բ���䵽����ں� fin��ֻ��Ϊ���򣱣���Լ�����ƣ�
        % ------------------------------------- %
        % ��sum_fij����ԭ����sum(fij{t-1}(:,j))
        sum_fij = 0;
        % ��sum_Fij����ԭ����sum(Fij{t-1}(:,j))
        sum_Fij = 0;
        for ind=1:size(conflict_fij{t-1}{j}, 1)           
            sum_fij = sum_fij + fij{t-1}( conflict_fij{t-1}{j}(ind,1), conflict_fij{t-1}{j}(ind,2) );
            sum_Fij = sum_Fij + Fij{t-1}( conflict_fij{t-1}{j}(ind,1), conflict_fij{t-1}{j}(ind,2) );
        end
        % ------------------------------------- %
        all_fin = sum_fij + fsj{t}(j) + sum(fmj{t}(j,:)) + sum_fid + sum_fiv;        
        % ����Բ�ı�׼�� ��ں� Fin����ͬ�ϣ�
        all_Fin = sum_Fij + Fsj{t}(j) + sum(Fmj{t}(j,:)) + sum_Fid + sum_Fiv;
        % ------------------------------------- %   
        
        if 1-all_Fin == 1 % ����GT���龰���ֵĴ���
            fd_Tcount = fd_Tcount + 1;
        end
        if 1-all_fin == 1 % ����������龰���ֵĴ���
            fd_Pcount = fd_Pcount + 1;
        end
        % ������ʧ������ֻ�е�����Բ��ʵ��ں�Ϊ0�����龰������������ں�Ϊ1ʱ����������ʧ
        fd_TP = fd_TP + (1 - all_Fin)*(1 - all_fin);
        
    end
end

% ===========���ǵ�һ֡�ĳ��ڣ�����Բj��ʵ����Ϊ0���ҷ������Ϊ1�������һ����ʧ
t = s_frame;
for j=1:n(t)
    sum_fmj = 0;
    sum_Fmj = 0;
    
    % sum_fmj Ϊ���а��� j ���ں� pair �� fmj ֮��
    for ind=1:numel(conflict_pair_last_xy{t}{j})/2
        sum_fmj = sum_fmj + fmj{t+1}( conflict_pair_last_xy{t}{j}(ind,1), conflict_pair_last_xy{t}{j}(ind,2) );
        sum_Fmj = sum_Fmj + Fmj{t+1}( conflict_pair_last_xy{t}{j}(ind,1), conflict_pair_last_xy{t}{j}(ind,2) );
    end
    % ����Բ���䵽�ĳ��ں� fout��ֻ��Ϊ���򣱣���Լ�����ƣ�
    all_fout = sum(fij{t}(j,:)) + fit{t}(j) + sum(fid{t}(j,:)) + sum(fiv{t}(j,:)) + sum_fmj;
    % ����Բ�ı�׼�𰸳��ں� Fout����ͬ�ϣ�
    all_Fout = sum(Fij{t}(j,:)) + Fit{t}(j) + sum(Fid{t}(j,:)) + sum(Fiv{t}(j,:)) + sum_Fmj;
            
    if 1-all_Fout == 1 % ����GT���龰���ֵĴ���
        fd_Tcount = fd_Tcount + 1;
    end
    if 1-all_fout == 1 % ����������龰���ֵĴ���
        fd_Pcount = fd_Pcount + 1;
    end

    % ������ʧ������ֻ�е�����Բ��ʵ���ں�Ϊ0�����龰������������ں�Ϊ1ʱ����������ʧ
    fd_TP = fd_TP + (1 - all_Fout)*(1 - all_fout);
end

% ��Ϊ�ж��˵�Ĵ��ڣ����龰�����ִ����ز�Ϊ0
% ע�⣬�˴���ָ�ġ��龰�������� �������龰 �� ��ì�ܼ��ų�������Բ��ȷ������Ϊ�����/�����Ƿ�Ϊ0
    
%% 5.�õ�Ŀ�꺯��
% ======================================================================= %
%
% �ܽ�Ŀ�꺯������ɣ�w'* fai(x,z) + delta(z,z*)
% ���� delta(z,z*) ����2������ɣ�
% (1) F*(1 - f)     �����еĶ���    1 0ʱ������ʧ ����ÿ���¼������̱���f
%     ��6��: TN.fij = 0;
%            TN.fit = 0;
%            TN.fid = 0;
%            TN.fiv = 0;
%            TN.fmj = 0;
%            TN.fsj = 0;
%
% (2) (1 - F)* f	���Լ��Ķ���   0 1ʱ������ʧ ֻ������¼����̱���ȫ0����Բ
%     ��2��: 2�����֡�����cost�͵�һ֡�ĳ���cost
%            ���� add_cost
%
% ======================================================================= %
% ȫ����ʧ����cost֮��
% ���Լ��ķ�������ע�͵���
% sum_cost = TN.fij + TN.fit + TN.fid + TN.fiv + TN.fmj + TN.fsj + add_cost;
% ���ķ���

event_FN = sum(Tcount) - event_TP; % ϸ���¼���FN������û������ʵ�ʷ����ˣ�
fd_FN = fd_Tcount - fd_TP; % �龰��FN�����Բ�Ϊ�鵫ʵ��Ϊ�飩
cost_for_train = (event_FN + fd_FN)/ sum(Tcount); % ����recall�ļ��㣨������SSVMѵ��ʱ��cost��

% ------------------------ �������ڼ�����Լ����� ------------------------- %
Preci.all = (event_TP+fd_TP)/(sum(Pcount)+ fd_Pcount);
Preci.move = TP.fij/Pcount(1); % precision=TP/(TP+FP);
Preci.disappear = TP.fit/Pcount(2);
Preci.divide = TP.fid/Pcount(3);
Preci.split = TP.fiv/Pcount(4);
Preci.merge = TP.fmj/Pcount(5);
Preci.appear = TP.fsj/Pcount(6);
Preci.false_detection = fd_TP/fd_Pcount;

Recall.all = (event_TP+fd_TP)/(sum(Tcount)+ fd_Tcount);
Recall.move = (TP.fij)/Tcount(1); % recall=TP/(TP+FN)
Recall.disappear = (TP.fit)/Tcount(2);
Recall.divide = (TP.fid)/Tcount(3);
Recall.split = (TP.fiv)/Tcount(4);
Recall.merge = (TP.fmj)/Tcount(5);
Recall.appear = (TP.fsj)/Tcount(6);
Recall.false_detection = fd_TP/fd_Tcount;

FM.all = 1/(1/Preci.all + 1/Recall.all)*2;
FM.move = 1/(1/Preci.move + 1/Recall.move)*2;
FM.disappear = 1/(1/Preci.disappear + 1/Recall.disappear)*2;
FM.divide = 1/(1/Preci.divide + 1/Recall.divide)*2;
FM.split = 1/(1/Preci.split + 1/Recall.split)*2;
FM.merge = 1/(1/Preci.merge + 1/Recall.merge)*2;
FM.appear = 1/(1/Preci.appear + 1/Recall.appear)*2;
FM.false_detection = 1/(1/Preci.false_detection + 1/Recall.false_detection)*2;

PRF.Preci = Preci; % ��PRF���������еľ��������
PRF.Recall = Recall;
PRF.FM = FM;

COUNT.Tcount = Tcount; % ��COUNT���������еļ��������
COUNT.Pcount = Pcount;
COUNT.fd_Tcount = fd_Tcount;
COUNT.fd_Pcount = fd_Pcount;

% ----------------------------------------------------------------------- %


























