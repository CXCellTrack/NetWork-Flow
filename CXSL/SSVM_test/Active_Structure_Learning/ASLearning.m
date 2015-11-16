% ================================================================== %
%
% CX 2015.10.5
% 这个函数用于实现AcitiveLearning中的2帧2帧测试的做法，由脚本 CX_AL_Results 调用
% 将每次2帧的分配结果返回，用于计算最终的精度等信息
%
% ================================================================== %

function [ fij fit fid fiv fmj fsj F ] = ASLearning( s_frame, e_frame )

% 指定测试帧的范围
% s_frame = 1;
% e_frame = 2;

% 指定在哪个数据集上进行计算（training or competition）
dataset = 'competition';
disp('  ============================');
disp(['  计算 ',num2str(s_frame), '—',num2str(e_frame), ' 帧的约束条件...']);
% 组建目标函数 注意：不包含损失函数
% ----------------------------------------------------------------------- %
% B方法：先计算 prob = <w,feature>（41s）在计算 obj = prob.*z（1s）
%        如果放在循环中，每次花一秒计算obj略长，但只计算一次的话速度很快
% 因此此处使用B方法速度较快，经验证 count_F_false 的计算没有问题

[ fij fit fid fiv fmj fsj ] = CXSL_Assign_FlowVar( dataset, s_frame, e_frame );
% 此处的true/false决定是否加入可选约束（要与训练时的选择一致）
[ F ] = CXSL_Calculate_Constraint_New_Conflict( dataset, [3 5], s_frame, e_frame, fij, fit, fid, fiv, fmj, fsj);
% 计算目标函数（需要载入之前计算好的特征）
% object_function = CXSL_Calculate_Obj_New( dataset, w_best, s_frame, e_frame, fij, fit, fid, fiv, fmj, fsj );

% ----------------------------------------------------------------------- %

%% 最终求解
% disp('  开始求解ILP...');
% % clearvars -except F object_function s_frame e_frame  fij fid fiv fit fsj fmj loss dataset count count_F_false exist_GT;
% % 注意，原先采用先算出 fai(x,z) = <feature,z>，在计算 obj = <w,fai(x,z)>;
% % 现在采用先计算 <w,feature>，在计算 obj = <w,feature>*z，速度得到了明显提升
% % 但这是针对于一次计算而言，如果在循环中每次都要这么计算目标函数，速度还是没有原方法快 2015.6.24
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








