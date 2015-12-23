% ======================================================================= %
%
% ����ǲ��� active structured learning �з������ȵĽű�
% ����2֡����������ѵ�����w��22���ԣ��õ����䷽��F���������frame-1��
% 2015.10.6
% ======================================================================= %

clear;close all

% ָ���Ƿ����GT��������ڣ�����㾫�ȵ�ָ�꣬������Ҫ��
exist_GT = 1;
dataset = 'competition';
[ ~, trackpath ] = getpath( dataset );
[ ~, traintrackpath ] = getpath( 'training' );

%% �����22֮���Լ��
global Ellipse n conflict_fij conflict_pair_next_xy conflict_pair_last_xy
load([ trackpath, '\Pair\Pre_data_New.mat']);
frame = numel(n);

fij = cell(frame-1,1);
fid = cell(frame-1,1);
fiv = cell(frame-1,1);
fit = cell(frame-1,1);
fsj = cell(frame,1);
fmj = cell(frame,1);
F = {};

% Ԥ����Լ������
for t=1:frame-1
    s_frame = t;
    e_frame = t + 1;
    % ʹ��xlou�� sl �ķ�������ѵ��
%     [ Fij, Fit, Fid, Fiv, Fsj, Fmj ] = ASLearning( w_best, s_frame, e_frame, Fij, Fit, Fid, Fiv, Fsj, Fmj );
    [ fij{t} fit{t} fid{t} fiv{t} fmj{t} fsj{t} F{t} ] = ASLearning( s_frame, e_frame );
end

%% ����w����Ŀ�꺯��������⣨�����ֿ������ڸ���w��
Fij = cell(frame-1,1);
Fid = cell(frame-1,1);
Fiv = cell(frame-1,1);
Fit = cell(frame-1,1);
Fsj = cell(frame,1);
Fmj = cell(frame,1);

if 0
    disp('  ����֮ǰ SSVM ѵ���õ������ w...');
    load([ traintrackpath, '\�ṹ��ѧϰ\SSVM_Best_W_New.mat']);
else % ��Ҫ���� excel �м��ص���ǰ��ʵ������ֻ��Ҫ�ֶ���д w_best ����
    disp('  �����ֶ���д��w...');
    % ע�⣡ֻ����2֡����������ѵ����������������  
    thisfile = 'BCFWavg_paper\40_2\loss_40_2_cons35_cost1_initwp_line.mat';
    load([ traintrackpath, '\ѵ�������¼\', thisfile ], 'w_best','Wavg');
%     w_best = Wavg{1792};
end

options = sdpsettings('verbose',0,'solver','gurobi');
for t=1:frame-1
    s_frame = t;
    e_frame = t + 1;
    object_function = CXSL_Calculate_Obj_New( dataset, w_best, s_frame, e_frame, fij{t}, fit{t}, fid{t}, fiv{t}, fmj{t}, fsj{t} );
    % ���
    disp(['  ��� ',num2str(s_frame), '��',num2str(e_frame), ' ֡��Ŀ�꺯��...']);
    sol = solvesdp( F{t}, -object_function, options );
    if sol.problem == 0
        for i = s_frame:e_frame-1
            Fij{i} = round(value(fij{t}{i})) ;
            Fid{i} = round(value(fid{t}{i})) ;
            Fiv{i} = round(value(fiv{t}{i})) ;
            Fit{i} = round(value(fit{t}{i})) ;
        end
        for i = s_frame+1:e_frame
            Fsj{i} = round(value(fsj{t}{i})) ;
            Fmj{i} = round(value(fmj{t}{i})) ;
        end

        COST = value(object_function);
        fprintf('\tcost:\t%.4f\n\n', COST);
        % ------------------------------------------------------ %
    else
        sol.info
        yalmiperror(sol.problem)
    end
end

%% ���ú������㾫�ȣ���˵���ȣ�
s_frame = 1;
e_frame = frame;
addfd = 0; % �����ܾ����м����龰
[ ~, PRF, COUNT ] = CX_Calculate_Loss( dataset, addfd, exist_GT, s_frame, e_frame, Fij, Fit, Fid, Fiv, Fmj, Fsj );

if isa(PRF, 'struct')
    
    disp(PRF.FM) % ֻ��ʾF-measure������   
    % ��ӡ���¼�����
    Preci = struct2array(PRF.Preci)*100;
    Recall = struct2array(PRF.Recall)*100;
    FMeasure = struct2array(PRF.FM)*100;
    
    PRM_for_excel = [ Preci;Recall;FMeasure ];
    COUNT_for_excel = [ COUNT.Tcount', COUNT.fd_Tcount; COUNT.Pcount', COUNT.fd_Pcount ];
    disp(COUNT_for_excel)
end
% ------------------------------------------------------ %

%% ����õ���������������� track_data �У���������ͼ�ã�ͨ����Ҫ���棡��
if 0
    if exist('thisfile', 'var')
        matpath = [trackpath, '\���Խ����¼\',strrep(thisfile, 'loss_', '')];
    else
        matpath = [trackpath, '\���Խ����¼\local.mat'];
    end
    save(matpath, 'PRF','COUNT','Fij','Fit','Fid','Fiv','Fmj','Fsj'); % ע���޸�mat����
    file = fopen(strrep(matpath,'mat','txt'), 'w'); fclose(file);
end










