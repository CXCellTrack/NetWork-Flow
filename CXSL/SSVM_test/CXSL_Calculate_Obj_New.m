function object_function = CXSL_Calculate_Obj_New( dataset, w_best, s_frame, e_frame, fij, fit, fid, fiv, fmj, fsj )

% ======================================================================= %
%
% �������������debug���亯��ԭ������Ϊ CXSL_ILP_Using_Best_W �ڲ�����
% ������ʵ����һ�ּ���Ŀ�꺯���ķ������� <w,feature>.*z
%
% ======================================================================= %

% ----------------------------------------------------------------------- %
[ ~, trackpath ] = getpath( dataset );
[ ~, traintrackpath ] = getpath( 'training' );
load([ traintrackpath, '\�ṹ��ѧϰ\initial_w_New.mat']);
df = zeros(6,1);
df(1) = 1;
df(2) = numel(wij)+1;
df(3) = numel(wij)+numel(wit)+1;
df(4) = numel(wij)+numel(wit)+numel(wid)+1;
df(5) = numel(wij)+numel(wit)+numel(wid)+numel(wiv)+1;
df(6) = numel(wij)+numel(wit)+numel(wid)+numel(wiv)+numel(wmj)+1;

wij = w_best( df(1):df(2)-1 )'; % ��̬�ָֵ���Ժ��޸������������޸�������
wit = w_best( df(2):df(3)-1 )';
wid = w_best( df(3):df(4)-1 )';
wiv = w_best( df(4):df(5)-1 )';
wmj = w_best( df(5):df(6)-1 )';
wsj = w_best( df(6):end )';

wij = w_best( 1:8 )'; % ���������w
wit = w_best( 9:13 )';
wid = w_best( 14:21 )';
wiv = w_best( 22:25 )';
wmj = w_best( 26:29 )';
wsj = w_best( 30:34 )';








% ��������
load([ trackpath, '\�ṹ��ѧϰ\Feature_New.mat']);
% ----------------------------------------------------------------------- %

obj_ij = 0;
obj_it = 0;
obj_id = 0;
obj_iv = 0;
obj_mj = 0;
obj_sj = 0;

for tt=s_frame:e_frame-1
    % <w,feature>�պõõ�һ��������fij��˼��ɣ��ٶȼӿ�
    mat_fij = cellfun(@(x)wij*x, feature_fij{tt});
    obj_ij = obj_ij + sum(sum( mat_fij.*fij{tt} ));
    
    mat_fit = cellfun(@(x)wit*x, feature_fit{tt});
    obj_it = obj_it + sum(sum( mat_fit.*fit{tt} ));
    
    mat_fid = cellfun(@(x)wid*x, feature_fid{tt});
    obj_id = obj_id + sum(sum( mat_fid.*fid{tt} ));
    
    mat_fiv = cellfun(@(x)wiv*x, feature_fiv{tt});
    obj_iv = obj_iv + sum(sum( mat_fiv.*fiv{tt} ));
    
end
for tt=s_frame+1:e_frame
    
    mat_fmj = cellfun(@(x)wmj*x, feature_fmj{tt});
    obj_mj = obj_mj + sum(sum( mat_fmj.*fmj{tt} )) ;
    
    mat_fsj = cellfun(@(x)wsj*x, feature_fsj{tt});
    obj_sj = obj_sj + sum(sum( mat_fsj.*fsj{tt} ));
    
end
    
object_function = obj_ij + obj_it + obj_id + obj_iv + obj_mj + obj_sj;

end



