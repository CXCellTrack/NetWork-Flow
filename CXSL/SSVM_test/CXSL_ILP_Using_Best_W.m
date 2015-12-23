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
disp(['  ���� ',num2str(s_frame), '��',num2str(e_frame), ' ֡��Ŀ�꺯����Լ������...']);

%% ����Լ������ ע�⣺��������ʧ����
% ----------------------------------------------------------------------- %
% B�������ȼ��� prob = <w,feature>��41s���ڼ��� obj = prob.*z��1s��
%        �������ѭ���У�ÿ�λ�һ�����obj�Գ�����ֻ����һ�εĻ��ٶȺܿ�
% ��˴˴�ʹ��B�����ٶȽϿ죬����֤ count_F_false �ļ���û������
disp('�������̱���...');tic
[ fij fit fid fiv fmj fsj ] = CXSL_Assign_FlowVar( dataset, s_frame, e_frame );
toc;disp('����Լ������...');
% �˴���true/false�����Ƿ�����ѡԼ����Ҫ��ѵ��ʱ��ѡ��һ��)
use_op_cons_test = [3 5 4];
[ Ffull, Fbase ] = CXSL_Calculate_Constraint_New_Conflict( dataset, use_op_cons_test, s_frame, e_frame, fij, fit, fid, fiv, fmj, fsj);
% ����Ŀ�꺯������Ҫ����֮ǰ����õ�������
if 1
    F = Ffull;
else
    F = Fbase; % ���ÿ�ѡԼ��
end

%% �齨Ŀ�꺯��
[ ~, traintrackpath ] = getpath( 'training' );
if 1
    if 0
        disp('  ����֮ǰ SSVM ѵ���õ������ w...');
        load([ traintrackpath, '\�ṹ��ѧϰ\SSVM_Best_W_New.mat']);
    else % ��Ҫ���� excel �м��ص���ǰ��ʵ������ֻ��Ҫ�ֶ���д w_best ����
        disp('  �����ֶ���д��w...');
        thisfile = 'BCFWavg_paper\withgap\64_2\loss_64_2_cons35_cost1_initwp_line_b_rng.mat';
        load([ traintrackpath, '\ѵ�������¼\', thisfile ], 'w_best','use_op_cons','Wavg');
%         w_best = Wavg{394};
%         if ~isequal(use_op_cons, use_op_cons_test)
%             error('�������õĿ�ѡԼ����ѵ����һ�£�');
%         end
    end
else  
    disp('  ����local SVM ѵ�������ĸ��¼�w...');
    % Ҳ����ʹ�õ���svmѵ�������ĸ��¼�w���ڽ������
    load([ traintrackpath, '\�ṹ��ѧϰ\initial_w_New.mat']);
    % ��ƫ��
    w_best = [ wij,bij, wit,bit, wid,bid, wiv,biv, wmj,bmj, wsj,bsj ]';
    % ��һ�����ݼ�
%     w_best = [-24.05027785	0.965959752	-0.209700235	0.023655542	-0.901444678	0.915527485	-0.723055368	0.78127324	-22.34659216	-3.45491283	-1.682414322	-5.355960441	-2.391659001	2.862181421	-7.382944338	8.382838223	1.94377663	-0.451290137	-1.07738777	-4.844423375	-1.122913059	-0.801496889	3.907101647	-11.61160994	3.710115534	0.998335816	4.252699702	0.790594494	1.207125853	3.799458373	1.390618031	5.18991389	1.129864864	0.673380786	-2.076937813	-1.97433464	-1.980221778	-0.051210814	0.597328997	-3.897482158]';
    % �ڶ������ݼ�
%     w_best = [-18.21315239	2.048551055	0.090611096	-0.07830978	-0.681768441	0.091705287	-0.284558766	0.113666465	-16.77170209	-2.820207584	-1.606735489	-2.170929556	-2.000511632	2.42450433	-4.406444861	10.34098417	1.758814312	-0.819906672	-2.159095585	-4.969572233	1.29607646	0.202318113	3.379177651	-14.25716614	7.109097539	3.366559674	2.10659084	-2.499899814	-3.9849466	2.501397849	1.351917771	4.699879458	-0.525793863	-0.261628668	-4.945586679	-1.846739083	-4.998199696	-0.003648138	1.737582541	-8.35761154]';
    % ���������ݼ�
%     w_best = [-17.69956928	1.61055833	-0.001787913	0.199592206	-0.287005286	-1.058810184	-0.194715268	0.170821175	-15.86943474	-2.459629861	-1.102503407	-1.233255838	-1.209682446	-2.849135318	-5.932152347	3.635597276	3.112043952	0.905772748	-0.112320949	-0.658275898	4.373948249	0.504445114	-0.211398811	-3.571706443	0	0	0	0	0	0	0	0	0	0	-2.873279295	-0.074622812	-1.438747225	-0.381464188	-2.838608426	-5.896985464]';
end

disp('�齨Ŀ�꺯��...');tic
object_function = CXSL_Calculate_Obj_New( dataset, w_best, s_frame, e_frame, fij, fit, fid, fiv, fmj, fsj );
toc

%% �������
disp('  ��ʼ���ILP...');
% clearvars -except F object_function s_frame e_frame  fij fid fiv fit fsj fmj loss dataset count count_F_false exist_GT trackpath thisfile use_op_cons_test;
% ע�⣬ԭ�Ȳ�������� fai(x,z) = <feature,z>���ڼ��� obj = <w,fai(x,z)>;
% ���ڲ����ȼ��� <w,feature>���ڼ��� obj = <w,feature>*z���ٶȵõ�����������
% �����������һ�μ�����ԣ������ѭ����ÿ�ζ�Ҫ��ô����Ŀ�꺯�����ٶȻ���û��ԭ������ 2015.6.24

options = sdpsettings('verbose',0,'solver','gurobi');
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
fprintf('\tcost:\t%.4f\n\n', COST); % save([trackpath, '\�ṹ��ѧϰ\Tracking_Data.mat'], 'Fij','Fid','Fiv','Fit','Fmj','Fsj');

%% ���ú������㾫�ȣ���˵���ȣ�
addfd = 0; % ѡ���Ƿ��龰�������ܾ�����
% ָ���Ƿ����GT��������ڣ�����㾫�ȵ�ָ�꣬������Ҫ��
exist_GT = 1;
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
        matpath = [trackpath, '\���Խ����¼\local_b.mat'];
    end
    save(matpath, 'PRF','COUNT','Fij','Fit','Fid','Fiv','Fmj','Fsj'); % ע���޸�mat����
    file = fopen(strrep(matpath,'mat','txt'), 'w'); fclose(file);
end




