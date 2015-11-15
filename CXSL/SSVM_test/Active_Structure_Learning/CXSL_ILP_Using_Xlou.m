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


% ���ú������㾫�ȣ���˵���ȣ�
s_frame = 1;
e_frame = frame;
[ ~, PRF, COUNT ] = CX_Calculate_Loss( dataset, exist_GT, s_frame, e_frame, Fij, Fit, Fid, Fiv, Fmj, Fsj );

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
    txtpath = [ trackpath, '\���Խ����¼\BCFWavg_my\40_2_cons1235.txt'];
    fidin = fopen(txtpath, 'w'); fclose(fidin);
    save(strrep(txtpath,'txt','mat'), 'PRF','COUNT','Fij','Fid','Fiv','Fit','Fsj','Fmj'); % ע���޸�mat����
    
%     save([ trackpath, '\�ṹ��ѧϰ\Tracking_Data.mat'], 'Fij','Fid','Fiv','Fit','Fsj','Fmj');
end










