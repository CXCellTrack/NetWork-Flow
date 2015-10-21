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

if 1
    if 1
        disp('  ����֮ǰ SSVM ѵ���õ������ w...');
        load([ traintrackpath, '\�ṹ��ѧϰ\SSVM_Best_W_New.mat']);
    else % ��Ҫ���� excel �м��ص���ǰ��ʵ������ֻ��Ҫ�ֶ���д w_best ����
        disp('  �����ֶ���д��w...');
        % ע�⣡ֻ����2֡����������ѵ����������������
        w_best = [-5.170345124	0.574560068	0.436206287	0.258141338	-0.172298816	0.080163188	0.024041306	0.089464995	-4.252191081	-0.274891168	-1.052420984	-0.330896797	-0.36331268	1.32905016	-1.735174617	3.023747434	1.17704628	0.011955574	-0.463793701	-1.254599365	0.162414663	0.181850867	1.291465386	-3.83609205	1.852794508	0.903066428	-0.199206694	-0.390306403	-1.321985405	0.608154159	0.347836237	1.289479594	0.221317699	-0.439501761	-1.135850144	-0.008892086	-1.223631805	-0.153622387	0.195378585	-2.347423748]';
    end
else  
    disp('  ���뵥�� SVM ѵ�������ĸ��¼�w���ٽ������...');
    % Ҳ����ʹ�õ���svmѵ�������ĸ��¼�w���ڽ������
    load([ traintrackpath, '\�ṹ��ѧϰ\initial_w_New.mat']);
    % ע�� w ��˳������
    w_best = [ wij, wit, wid, wiv, wmj, wsj ]';
end

load([ trackpath, '\Pair\Pre_data_New.mat'], 'n');
frame = numel(n);
Fij = cell(frame-1,1);
Fid = cell(frame-1,1);
Fiv = cell(frame-1,1);
Fit = cell(frame-1,1);
Fsj = cell(frame,1);
Fmj = cell(frame,1);

% ����cplex����22���
for i=1:frame-1
    s_frame = i;
    e_frame = i + 1;
    % ʹ��xlou�� sl �ķ�������ѵ��
    [ Fij, Fit, Fid, Fiv, Fsj, Fmj ] = ASLearning( w_best, s_frame, e_frame, Fij, Fit, Fid, Fiv, Fsj, Fmj );
end

fij = Fij;
fid = Fid;
fiv = Fiv;
fit = Fit;
fsj = Fsj;
fmj = Fmj;
% ���ú������㾫�ȣ���˵���ȣ�
s_frame = 1;
e_frame = frame;
[ ~, PRF, COUNT ] = CX_Calculate_Loss( dataset, exist_GT, s_frame, e_frame, fij, fit, fid, fiv, fmj, fsj );

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


% ����õ���������������� track_data �У���������ͼ�ã�ͨ����Ҫ���棡��
if 0
    txtpath = [ trackpath, '\���Խ����¼\BCFWavg_New\32_2_y.txt'];
    fid = fopen(txtpath, 'w'); fclose(fid);
    save(strrep(txtpath,'txt','mat'), 'PRF','COUNT','fij','fid','fiv','fit','fsj','fmj'); % ע���޸�mat����
    
%     save([ trackpath, '\�ṹ��ѧϰ\Tracking_Data.mat'], 'Fij','Fid','Fiv','Fit','Fsj','Fmj');
end










