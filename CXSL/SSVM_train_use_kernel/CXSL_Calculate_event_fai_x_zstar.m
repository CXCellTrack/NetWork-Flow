function [ phi_x_z ] = CXSL_Calculate_event_fai_x_zstar( s_frame, e_frame, hat_or_star)

% ============================= %
% ����������ڼ��� fai(x,z*) 
% Ҳ���Լ��� fai(x,z^) �� fai(x,z^) ������ILP�еĽ�� value ��������������ظ�����
% ============================= %
% ��������
dataset = 'training';
[ ~, trackpath ] = getpath( dataset );
load([ trackpath, '\�ṹ��ѧϰ\Feature_Plus_New.mat']);
% �����׼��GT
load([ trackpath, '\GT\GT_after_hand_tune\GT_Flow_Variables_New.mat']);

fai_fij = 0;
fai_fit = 0;
fai_fid = 0;
fai_fiv = 0;
fai_fmj = 0;
fai_fsj = 0;

if strcmp(hat_or_star, 'hat') % �ⲿ���Ѿ���ILP�Ľ����value�����������¼���
    % ���㵱ǰw�µ�fai(x,z^)
    for t = s_frame:e_frame-1   
        fai_fij = fai_fij + cell_dot_mutil( feature_fij_p{t}, Fij_c{t} );
        fai_fit = fai_fit + cell_dot_mutil( feature_fit_p{t}, Fit_c{t} );
        fai_fid = fai_fid + cell_dot_mutil( feature_fid_p{t}, Fid_c{t} );
        fai_fiv = fai_fiv + cell_dot_mutil( feature_fiv_p{t}, Fiv_c{t} );   
    end
    for t = s_frame+1:e_frame
        fai_fmj = fai_fmj + cell_dot_mutil( feature_fmj_p{t}, Fmj_c{t} );  
        fai_fsj = fai_fsj + cell_dot_mutil( feature_fsj_p{t}, Fsj_c{t} );
    end
    
elseif strcmp(hat_or_star, 'star')
    % �����׼���е�fai(x,z*)
%     tic;
    for t = s_frame:e_frame-1  
        
        fai_fij = fai_fij + cell_dot_mutil( feature_fij_p{t}, Fij{t} );
            
        fai_fit = fai_fit + cell_dot_mutil( feature_fit_p{t}, Fit{t} ); %  �Ӳ���sum = ��һ��
        
        fai_fid = fai_fid + cell_dot_mutil( feature_fid_p{t}, Fid{t} );  
		
        fai_fiv = fai_fiv + cell_dot_mutil( feature_fiv_p{t}, Fiv{t} );         
    end
    
    for t = s_frame+1:e_frame
	
        fai_fmj = fai_fmj + cell_dot_mutil( feature_fmj_p{t}, Fmj{t} );
		
        fai_fsj = fai_fsj + cell_dot_mutil( feature_fsj_p{t}, Fsj{t} );
    end
%     toc;
end

% �����е� fai_f ��ϳ�������������������
% fai_x_z = [ fai_fij; fai_fit; fai_fid; fai_fiv; fai_fmj; fai_fsj; ];
% ��ʹ�ú�ʱ�����÷ֿ��¼��� phi �Ƚ����
phi_x_z = cell(6,1);
phi_x_z{1} = fai_fij;
phi_x_z{2} = fai_fit;
phi_x_z{3} = fai_fid;
phi_x_z{4} = fai_fiv;
phi_x_z{5} = fai_fmj;
phi_x_z{6} = fai_fsj;

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
















