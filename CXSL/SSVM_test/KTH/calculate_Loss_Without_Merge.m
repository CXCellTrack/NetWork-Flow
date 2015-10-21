function [ PRF COUNT ] = calculate_Loss_Without_Merge( dataset, exist_GT, s_frame, e_frame, Pcount, fij, fit, fid, fsj )

% ���㾫��
%% ��������

[ ~, trackpath ] = getpath( dataset );

load([ trackpath, '\Pair\Pre_data_New.mat'], 'n','conflict_pair_last_xy','conflict_pair_next_xy','conflict_fij');
gt_flow_addr = [ trackpath, '\GT\GT_after_hand_tune\GT_Flow_Variables_New.mat'];

% �����׼��GT������Ҫ��GT���ڵ�����£�ʹ��bool exist_GT����ʶGT�Ƿ���ڣ�
if exist_GT
    load( gt_flow_addr );
else
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
TP.fsj = 0;

% PN������Խ����TF����GT���
%    T  F
% P  TP FP
% N  FN TN
% 
% precision = tp/(tp+fp)
% recall = tp/(tp+fn);

Tcount = zeros(4,1); % tp+tn
% ��ʧ���������� false negative ���㣬�� f*=1 && f=0 ʱ������ʧ
for t = s_frame:e_frame-1
    TP.fij = TP.fij + sum(sum( Fij{t}.*fij{t} ));
    Tcount(1) = Tcount(1) + sum(sum( Fij{t})); % ͳ�Ƴ� Fij{t}��1�ĸ���
%     Pcount(1) = Pcount(1) + sum(sum( fij{t})); % �����г��ֵĸ���
    
    TP.fit = TP.fit + sum(sum( Fit{t}.*fit{t} ));
    Tcount(2) = Tcount(2) + sum(sum( Fit{t}));
    
    TP.fid = TP.fid + sum(sum( Fid{t}.*fid{t} ));
    Tcount(3) = Tcount(3) + sum(sum( Fid{t}));
end
for t = s_frame+1:e_frame
    TP.fsj = TP.fsj + sum(sum( Fsj{t}.*fsj{t} ));
    Tcount(4) = Tcount(4) + sum(sum( Fsj{t}));
end

% ============= ע�⣺�����еķ����ǰ����е�count�������ڼ��� ============== %
event_TP = TP.fij + TP.fit + TP.fid + TP.fsj;
 
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

% ------------------------ �������ڼ�����Լ����� ------------------------- %
Preci.all = event_TP/(sum(Pcount));
Preci.move = TP.fij/Pcount(1); % precision=TP/(TP+FP);
Preci.disappear = TP.fit/Pcount(2);
Preci.divide = TP.fid/Pcount(3);
Preci.appear = TP.fsj/Pcount(4);

Recall.all = event_TP/(sum(Tcount));
Recall.move = (TP.fij)/Tcount(1); % recall=TP/(TP+FN)
Recall.disappear = (TP.fit)/Tcount(2);
Recall.divide = (TP.fid)/Tcount(3);
Recall.appear = (TP.fsj)/Tcount(4);

FM.all = 1/(1/Preci.all + 1/Recall.all)*2;
FM.move = 1/(1/Preci.move + 1/Recall.move)*2;
FM.disappear = 1/(1/Preci.disappear + 1/Recall.disappear)*2;
FM.divide = 1/(1/Preci.divide + 1/Recall.divide)*2;
FM.appear = 1/(1/Preci.appear + 1/Recall.appear)*2;

PRF.Preci = Preci; % ��PRF���������еľ��������
PRF.Recall = Recall;
PRF.FM = FM;

COUNT.Tcount = Tcount; % ��COUNT���������еļ��������
COUNT.Pcount = Pcount;

% ----------------------------------------------------------------------- %


























