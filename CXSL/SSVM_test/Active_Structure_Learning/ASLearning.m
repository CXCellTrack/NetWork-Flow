% ================================================================== %
%
% CX 2015.10.5
% �����������ʵ��AcitiveLearning�е�2֡2֡���Ե��������ɽű� CX_AL_Results ����
% ��ÿ��2֡�ķ��������أ����ڼ������յľ��ȵ���Ϣ
%
% ================================================================== %

function [ fij fit fid fiv fmj fsj F ] = ASLearning( s_frame, e_frame )

% ָ������֡�ķ�Χ
% s_frame = 1;
% e_frame = 2;

% ָ�����ĸ����ݼ��Ͻ��м��㣨training or competition��
dataset = 'competition';
disp('  ============================');
disp(['  ���� ',num2str(s_frame), '��',num2str(e_frame), ' ֡��Լ������...']);
% �齨Ŀ�꺯�� ע�⣺��������ʧ����
% ----------------------------------------------------------------------- %
% B�������ȼ��� prob = <w,feature>��41s���ڼ��� obj = prob.*z��1s��
%        �������ѭ���У�ÿ�λ�һ�����obj�Գ�����ֻ����һ�εĻ��ٶȺܿ�
% ��˴˴�ʹ��B�����ٶȽϿ죬����֤ count_F_false �ļ���û������

[ fij fit fid fiv fmj fsj ] = CXSL_Assign_FlowVar( dataset, s_frame, e_frame );
% �˴���true/false�����Ƿ�����ѡԼ����Ҫ��ѵ��ʱ��ѡ��һ�£�
[ F ] = CXSL_Calculate_Constraint_New_Conflict( dataset, [3 5], s_frame, e_frame, fij, fit, fid, fiv, fmj, fsj);
% ����Ŀ�꺯������Ҫ����֮ǰ����õ�������
% object_function = CXSL_Calculate_Obj_New( dataset, w_best, s_frame, e_frame, fij, fit, fid, fiv, fmj, fsj );

% ----------------------------------------------------------------------- %

%% �������
% disp('  ��ʼ���ILP...');
% % clearvars -except F object_function s_frame e_frame  fij fid fiv fit fsj fmj loss dataset count count_F_false exist_GT;
% % ע�⣬ԭ�Ȳ�������� fai(x,z) = <feature,z>���ڼ��� obj = <w,fai(x,z)>;
% % ���ڲ����ȼ��� <w,feature>���ڼ��� obj = <w,feature>*z���ٶȵõ�����������
% % �����������һ�μ�����ԣ������ѭ����ÿ�ζ�Ҫ��ô����Ŀ�꺯�����ٶȻ���û��ԭ������ 2015.6.24
% 
% options = sdpsettings('verbose',0,'solver','gurobi');
% sol = solvesdp( F, -object_function, options );
% 
% if sol.problem == 0
%     for t = s_frame:e_frame-1
%         Fij{t} = round(value(fij{t})) ;
%         Fid{t} = round(value(fid{t})) ;
%         Fiv{t} = round(value(fiv{t})) ;
%         Fit{t} = round(value(fit{t})) ;
%     end
%     for t = s_frame+1:e_frame
%         Fsj{t} = round(value(fsj{t})) ;
%         Fmj{t} = round(value(fmj{t})) ;
%     end
% 
%     COST = value(object_function);
%     fprintf('\tcost:\t%.4f\n\n', COST);
%     % ------------------------------------------------------ %
% else
%     sol.info
%     yalmiperror(sol.problem)
% end








