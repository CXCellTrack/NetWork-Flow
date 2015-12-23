function object_function = CXSL_Calculate_Obj_New( dataset, w_best, s_frame, e_frame, fij, fit, fid, fiv, fmj, fsj )

% ======================================================================= %
%
% �������������debug���亯��ԭ������Ϊ CXSL_ILP_Using_Best_W �ڲ�����
% ������ʵ����һ�ּ���Ŀ�꺯���ķ������� <w,feature>.*z
%
% ======================================================================= %

% ----------------------------------------------------------------------- %
[ ~, trackpath ] = getpath( dataset );
% ������ w_best Ϊ���������˴���ת��Ϊ�������ٽ������
% ��������
if numel(w_best)>=40
    load([ trackpath, '\�ṹ��ѧϰ\Feature_Plus_New.mat']);
else
    load([ trackpath, '\�ṹ��ѧϰ\Feature_New.mat']);
end

if exist('feature_fij_p','var')
    % ��Ϊ�����������ʹ�������w 
    d_f = [ 13, 12, 25, 24, 24, 12 ]; % ���¼�������ά��(Ҫ��˳��)
    d_f = [ 9, 6, 9, 5, 5, 6 ];
    tmp = cumsum(d_f);
    wij = w_best( 1:tmp(1) )'; 
    wit = w_best( tmp(1)+1:tmp(2) )';
    wid = w_best( tmp(2)+1:tmp(3) )';
    wiv = w_best( tmp(3)+1:tmp(4) )';
    wmj = w_best( tmp(4)+1:tmp(5) )';
    wsj = w_best( tmp(5)+1:tmp(6) )'; 

    feature_fij = feature_fij_p;
    feature_fit = feature_fit_p;
    feature_fid = feature_fid_p;
    feature_fiv = feature_fiv_p;
    feature_fmj = feature_fmj_p;
    feature_fsj = feature_fsj_p;
else
    % ���������w
    wij = w_best( 1:8 )'; 
    wit = w_best( 9:13 )';
    wid = w_best( 14:21 )';
    wiv = w_best( 22:25 )';
    wmj = w_best( 26:29 )';
    wsj = w_best( 30:34 )';
   
end
    
% ���Ŀ�꺯��
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
    % ������Щ����Ϊnan����������Բ��˵�������⣡������������ֲ���
    mat_fid(isnan(mat_fid)) = 0;
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



