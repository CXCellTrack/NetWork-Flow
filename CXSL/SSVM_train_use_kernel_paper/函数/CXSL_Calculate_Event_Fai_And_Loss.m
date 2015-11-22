function [ fij fit fid fiv fmj fsj phi_x_z cost_for_train cost_for_train_all ] = CXSL_Calculate_Event_Fai_And_Loss( s_frame, e_frame )

%
% �ڸ��� w ������£����ɱ��������Ŀ�꺯��
% ���ù�ʽ��2)��L(x,z,w) = w'*fai(x,z)
% ���� w' = [ wij, wit, wid, wiv, wmj, wsj ]  
% fai(x,z)Ҳ����ͬ��˳������
% �����������¼��ı��������Ŀ�꺯�� object_f���͵��ⲿ����ILP���
%
% s_frame Ϊ��ʼ֡  e_frame Ϊ����֡
% tic;
dataset = 'training'; % �������Ҳֻ��ѵ����ʹ��
[ ~, trackpath ] = getpath( dataset );

load([ trackpath, '\Pair\Pre_data_New.mat'], 'n','conflict_pair_last_xy','conflict_pair_next_xy','conflict_fij');
% ��������
load([ trackpath, '\�ṹ��ѧϰ\Feature_Plus_New.mat']);
% �����׼��GT
load([ trackpath, '\GT\GT_after_hand_tune\GT_Flow_Variables_New.mat']);

% ����1�ѱ�ע��
%% 2.������������
fij = cell(e_frame-1,1);
fid = cell(e_frame-1,1);
fit = cell(e_frame-1,1);
fiv = cell(e_frame-1,1);
fsj = cell(e_frame,1);
fmj = cell(e_frame,1);

for t = s_frame:e_frame-1
    %  t�е�m��ǰ���еĵ�width��ϸ��   
    fij{t} = binvar(n(t), 4, 'full'); %%fij��������
    fit{t} = binvar(n(t), 1, 'full');   %%��ʧ
    fid{t} = binvar(n(t), 6, 'full');   %%ĸϸ��
    fiv{t} = binvar(n(t), 6, 'full');   %%����
end

for t = s_frame+1:e_frame
    fmj{t} = binvar(n(t), 6,  'full');   %%�ں� ##ע������������t * t-1�����������ǣ�t * t+1��
    fsj{t} = binvar(n(t), 1, 'full');   %%���� 
end

%% 3.����ԭʼĿ�꺯��
% ���ù�ʽ��2)��L(x,z,w) = w'*fai(x,z)
% ���ù�ʽ��8): delta(z,z*) = 1/|z*| * sum(z* * (1-z))

fai_fij = 0;
fai_fit = 0;
fai_fid = 0;
fai_fiv = 0;
fai_fmj = 0;
fai_fsj = 0;

% ======================================================================= %
% ��������ȫ�������Ҫʹ�ò���1�еõ�ӳ������
% ��������ȫ�������ֱ�Ӿ����˸��죨6.5�� VS 31.4�룩
% ======================================================================= %
% A.ʹ�þ���˷����㣺
% tic;
for t = s_frame:e_frame-1
    % ȥ��sum��֮������ٶȸ��� �� 2.9168s �� 2.82s
    fai_fij = fai_fij + cell_dot_mutil( feature_fij_p{t}, fij{t} );

    fai_fit = fai_fit + cell_dot_mutil( feature_fit_p{t}, fit{t} );

    fai_fid = fai_fid + cell_dot_mutil( feature_fid_p{t}, fid{t} );

    fai_fiv = fai_fiv + cell_dot_mutil( feature_fiv_p{t}, fiv{t} );
end
% toc;
for t = s_frame+1:e_frame
    
    fai_fmj = fai_fmj + cell_dot_mutil( feature_fmj_p{t}, fmj{t} );

    fai_fsj = fai_fsj + cell_dot_mutil( feature_fsj_p{t}, fsj{t} );
end

% �����е� fai_f ��ϳ�������������������
% fai_x_z = [ fai_fij; fai_fit; fai_fid; fai_fiv; fai_fmj; fai_fsj; ];
% ��ʹ�ú�ʱ�����÷ֿ��¼��� phi �Ƚ����
phi_x_z = cell(6,1);
phi_x_z{1} = fai_fij;
phi_x_z{2} = fai_fit;
phi_x_z{3} = fai_fid;
phi_x_z{4} = fai_fiv;
phi_x_z{5} = fai_fmj;
phi_x_z{6} = fai_fsj;

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
    
    TP.fit = TP.fit + sum(sum( Fit{t}.*fit{t} ));
    Tcount(2) = Tcount(2) + sum(sum( Fit{t}));
    
    TP.fid = TP.fid + sum(sum( Fid{t}.*fid{t} ));
    Tcount(3) = Tcount(3) + sum(sum( Fid{t}));
    
    TP.fiv = TP.fiv + sum(sum( Fiv{t}.*fiv{t} ));
    Tcount(4) = Tcount(4) + sum(sum( Fiv{t}));
end
for t = s_frame+1:e_frame
    TP.fmj = TP.fmj + sum(sum( Fmj{t}.*fmj{t} ));
    Tcount(5) = Tcount(5) + sum(sum( Fmj{t}));
    
    TP.fsj = TP.fsj + sum(sum( Fsj{t}.*fsj{t} ));
    Tcount(6) = Tcount(6) + sum(sum( Fsj{t}));
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
cost_for_train_all = (event_FN + fd_FN)/ sum(Tcount); % ����recall�ļ��㣨������SSVMѵ��ʱ��cost��
cost_for_train = event_FN/ sum(Tcount);

end

function sum_f = cell_dot_mutil( feature, z )
%
% f Ϊ��������cell	n(t)*n(t+1)
% g Ϊ���̱�������
% ��� sum_f Ϊ��cell��������������֮��
% ===================================

% 1. for ѭ��ʽ���� 4.75��
% sum_f = zeros( size(feature{1,1}) );
% [h w] = size(feature);
% for i=1:h
%     for j=1:w
%         sum_f = sum_f + feature{i, j}* z(i,j);
%     end
% end

% 2. arrayfun���� 3.10��
sum_f = 0;
ss = numel(feature);
% feature1 = cellfun(@(x)x', feature, 'un',0);
feature2 = reshape(feature, ss, 1);
z1 = reshape(z, ss, 1);
f_z = arrayfun(@(x) feature2{x}*z1(x), 1:ss, 'un',0);
for x=1:ss
    sum_f = sum_f + f_z{x};
end

end











