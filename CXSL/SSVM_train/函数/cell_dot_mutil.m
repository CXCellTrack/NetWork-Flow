function sum_f = cell_dot_mutil( feature, z )
%
% ======================================================================= %
%
% ע�⣬�������ֻ����debug��ʵ���������е��õ���д�� CXSL_Calculate_Obj_With_Loss
% �� CXSL_Calculate_fai_x_zstar �е��ڲ�����
% �����Ѿ���֤�Ľ�����֣�
%       CXSL_Calculate_Obj_With_Loss ��������������Ϊsdpvar��ʹ�÷���������
%   ��arrayfun����
%       CXSL_Calculate_fai_x_zstar ����������Ϊdouble���ͣ�ʹ��forѭ��������
% ======================================================================= %

% feature Ϊ��������cell	n(t)*n(t+1)
% z Ϊ���̱�������
% ��� sum_f Ϊ��cell��������������֮��
% ===================================

% 1. for ѭ��ʽ���� 4.75��

% sum_f = zeros( size(feature{1,1}) );
% [h w] = size(feature);
% for i=1:h
%     for j=1:w
%         sum_f = sum_f + feature{i, j}* z(i,j);
%     end
% end

% 2. arrayfun���� 3.10��

% ������2���汾��һ���� CXSL_Calculate_Obj_With_Loss �м��� fai(x,z^) �İ汾
% ���ڵõ�����������Ϊ sdpvar���޷�����cell2mat��������һ��ֻ����sum���
% �������£�
% ================================================= %

sum_f = 0;
ss = numel(feature);
% feature1 = cellfun(@(x)x', feature, 'un',0);
feature2 = reshape(feature, ss, 1);
z1 = reshape(z, ss, 1);
f_z = arrayfun(@(x)feature2{x}*z1(x), 1:ss, 'un',0);
for x=1:ss
    sum_f = sum_f + f_z{x};
end

% ================================================= % 
% ��һ���汾�� CXSL_Calculate_fai_x_zstar �м��� fai(x,z*) �İ汾
% ������������Ϊdouble����˿���ת��Ϊmat��ʹ��sum��ͣ���������
% ��ʵ�ʹ����з������������û����һ���죬��˲�����

% ss = numel(feature);
% feature1 = cellfun(@(x)x', feature, 'un',0);
% feature2 = reshape(feature1, ss, 1);
% z1 = reshape(z, ss, 1);
% f_z = arrayfun(@(x)feature2{x}*z1(x), 1:ss, 'un',0)';
% f_z_mat = cell2mat(f_z);
% sum_f = sum(f_z_mat);

end



