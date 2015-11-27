% ================================================================== %
%
% CX 2015.6.23
% ����ű������ڸ��� w �������Ѱ�Ҳ��Լ�����ILP��������Ž⣬����ѷ��䷽��
% ��Ҫ�Լ�ָ����ʼ֡ s_frame �ͽ���֡ e_frame 
%
% ================================================================== %

clear;close all;

%% ָ�����ĸ����ݼ��Ͻ��м��㣨training or competition��
if 1
    dataset = 'competition';
    disp('<ѡ�����Լ�>');
else
    dataset = 'training';
    disp('<ѡ��ѵ����>');
end
[ segpath, trackpath ] = getpath( dataset );

% ָ������֡�ķ�Χ
segdir = dir([ segpath, '\*.tif']);
s_frame = 1;
e_frame = numel(segdir);
% ָ���Ƿ����GT��������ڣ�����㾫�ȵ�ָ�꣬������Ҫ��
exist_GT = 1;
disp(['  ���� ',num2str(s_frame), '��',num2str(e_frame), ' ֡��Ŀ�꺯����Լ������...']);

%% ����Լ������ ע�⣺��������ʧ����
% ----------------------------------------------------------------------- %
% B�������ȼ��� prob = <w,feature>��41s���ڼ��� obj = prob.*z��1s��
%        �������ѭ���У�ÿ�λ�һ�����obj�Գ�����ֻ����һ�εĻ��ٶȺܿ�
% ��˴˴�ʹ��B�����ٶȽϿ죬����֤ count_F_false �ļ���û������
disp('�������̱���...');tic
[ fij fit fid fiv fmj fsj ] = CXSL_Assign_FlowVar( dataset, s_frame, e_frame );
toc;disp('����Լ������...');
% �˴���true/false�����Ƿ�����ѡԼ����Ҫ��ѵ��ʱ��ѡ��һ�£�
[ F ] = CXSL_Calculate_Constraint_New_Conflict( dataset, true, s_frame, e_frame, fij, fit, fid, fiv, fmj, fsj);
% ����Ŀ�꺯������Ҫ����֮ǰ����õ�������

%% �齨Ŀ�꺯��
[ ~, traintrackpath ] = getpath( 'training' );
if 0
    disp('  ������һ�� SSVM ѵ���õ���w(���Ժ�) & A��aplha�������Ժˣ�...');
    load([ traintrackpath, '\�ṹ��ѧϰ\SSVM_Best_W_New.mat']);
else
    disp('  ����֮ǰ SSVM ѵ�������A��alpha...'); 
    name = 'loss_1_65_y';
    lossdir = [ traintrackpath, '\ѵ�������¼\BCFW\'];
    load([lossdir, name, '.mat']);
end

disp('�齨Ŀ�꺯��...');tic
object_function = CXSL_Calculate_Obj_hunhe( dataset, s_frame, e_frame,...
    w_best, A_best, alpha_best, kernel_type, cmd, islinear, ...
    fij, fit, fid, fiv, fmj, fsj );
toc

%% �������
disp('  ��ʼ���ILP...');
clearvars -except F object_function s_frame e_frame  fij fid fiv fit fsj fmj loss dataset count count_F_false exist_GT trackpath;
% ע�⣬ԭ�Ȳ�������� fai(x,z) = <feature,z>���ڼ��� obj = <w,fai(x,z)>;
% ���ڲ����ȼ��� <w,feature>���ڼ��� obj = <w,feature>*z���ٶȵõ�����������
% �����������һ�μ�����ԣ������ѭ����ÿ�ζ�Ҫ��ô����Ŀ�꺯�����ٶȻ���û��ԭ������ 2015.6.24
% options = sdpsettings('verbose',0,'solver','gurobi');
options = sdpsettings('verbose',0,'solver','cplex','saveduals',0);
sol = solvesdp( F, -object_function, options )

Fij = cell(e_frame-1,1);
Fid = cell(e_frame-1,1);
Fiv = cell(e_frame-1,1);
Fit = cell(e_frame-1,1);
Fsj = cell(e_frame,1);
Fmj = cell(e_frame,1);
if sol.problem ~= 0
    sol.info
    yalmiperror(sol.problem);
end

for t = s_frame:e_frame-1
    Fij{t} = round(value(fij{t})) ;
    Fid{t} = round(value(fid{t})) ;
    Fiv{t} = round(value(fiv{t})) ;
    Fit{t} = round(value(fit{t})) ;
end
for t = s_frame+1:e_frame
    Fsj{t} = round(value(fsj{t})) ;
    Fmj{t} = round(value(fmj{t})) ;
end

COST = value(object_function);
fprintf('\tcost:\t%.4f\n\n', COST);

% ���ú������㾫�ȣ���˵���ȣ�
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
    txtpath = [trackpath, '\���Խ����¼\BCFW\1_65_y.txt'];
    file = fopen(txtpath, 'w'); fclose(file);
    save(strrep(txtpath,'txt','mat'), 'PRF','COUNT','Fij','Fit','Fid','Fiv','Fmj','Fsj'); % ע���޸�mat����
    
%     save([ trackpath, '\�ṹ��ѧϰ\Tracking_Data.mat'], 'Fij','Fid','Fiv','Fit','Fsj','Fmj');
end





