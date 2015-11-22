function [ fai_x_z ] = CXSL_Calculate_phi_x_zstar_New( w, s_frame, e_frame, str)

% ============================= %
% ����������ڼ��� fai(x,z*) 
% Ҳ���Լ��� fai(x,z^) �� fai(x,z^) ������ILP�еĽ�� value ��������������ظ�����
% ============================= %

if isempty(w) % ��û������w��Ĭ��Ϊ��b��0
    w = zeros(40,1);
end

dataset = 'training';
[ ~, trackpath ] = getpath( dataset );
% ��������
if numel(w)==40
    load([ trackpath, '\�ṹ��ѧϰ\Feature_Plus_New.mat']);
    feature_fij = feature_fij_p;
    feature_fit = feature_fit_p;
    feature_fid = feature_fid_p;
    feature_fiv = feature_fiv_p;
    feature_fmj = feature_fmj_p;
    feature_fsj = feature_fsj_p;
else
    load([ trackpath, '\�ṹ��ѧϰ\Feature_New.mat']);
end
% �����׼��GT
load([ trackpath, '\GT\GT_after_hand_tune\GT_Flow_Variables_New.mat']);

fai_fij = 0;
fai_fit = 0;
fai_fid = 0;
fai_fiv = 0;
fai_fmj = 0;
fai_fsj = 0;

if strcmp(str, 'hat') % �ⲿ���Ѿ���ILP�Ľ����value�����������¼���
    %
end
if strcmp(str, 'star')
    % �����׼���е�fai(x,z*)
%     tic;
    for t = s_frame:e_frame-1  
        
        fai_fij = fai_fij + cell_dot_mutil( feature_fij{t}, Fij{t} );
            
        fai_fit = fai_fit + cell_dot_mutil( feature_fit{t}, Fit{t} ); %  �Ӳ���sum = ��һ��
        
        fai_fid = fai_fid + cell_dot_mutil( feature_fid{t}, Fid{t} );  
		
        fai_fiv = fai_fiv + cell_dot_mutil( feature_fiv{t}, Fiv{t} );         
    end
    
    for t = s_frame+1:e_frame
	
        fai_fmj = fai_fmj + cell_dot_mutil( feature_fmj{t}, Fmj{t} );
		
        fai_fsj = fai_fsj + cell_dot_mutil( feature_fsj{t}, Fsj{t} );
    end
%     toc;
end

% �����е� fai_f ��ϳ�������������������һ��42ά
fai_x_z = [ fai_fij; fai_fit; fai_fid; fai_fiv; fai_fmj; fai_fsj ];

end


%% ��������� CXSL_Calculate_Obj �еĺ���һ��
% ע�⣺��������Ϊdouble���ͣ�ʹ��forѭ��������

function sum_f = cell_dot_mutil( feature, z )
%
% f Ϊ��������cell	n(t)*n(t+1)
% g Ϊ���̱�������
% ��� sum_f Ϊ��cell��������������֮��
% ===================================

% 1. for ѭ��ʽ���� 0.035 ��

sum_f = zeros( size(feature{1,1}) );
[h w] = size(feature);
for i=1:h
    for j=1:w
        sum_f = sum_f + feature{i, j}* z(i,j);
    end
end

% 2. arrayfun���� 1.77 ��

% sum_f = 0;
% ss = numel(feature);
% % feature1 = cellfun(@(x)x', feature, 'un',0);
% feature2 = reshape(feature, ss, 1);
% z1 = reshape(z, ss, 1);
% f_z = arrayfun(@(x)feature2{x}*z1(x), 1:ss, 'un',0);
% for x=1:ss
%     sum_f = sum_f + f_z{x};
% end

end
















