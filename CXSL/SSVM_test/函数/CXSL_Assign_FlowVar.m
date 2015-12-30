function [ fij fit fid fiv fmj fsj ] = CXSL_Assign_FlowVar( dataset, s_frame, e_frame )

% ======================================================================= %
% ����������ڷ�����¼����̱�������������ʧ���� sum_cost
% �� CXSL_ILP_Using_Best_W ֱ�ӵ���������

% �ڸ��� w ������£����ɱ��������Ŀ�꺯��
% ���ù�ʽ��2)��L(x,z,w) = w'*fai(x,z)
% ���� w' = [ wij, wit, wid, wiv, wmj, wsj ]  
% fai(x,z)Ҳ����ͬ��˳������
% �����������¼��ı��������Ŀ�꺯�� object_f���͵��ⲿ����ILP���
%
% s_frame Ϊ��ʼ֡  e_frame Ϊ����֡

% ======================================================================= %
% tic;
[ ~, trackpath ] = getpath( dataset );
load([ trackpath, '\Pair\Pre_data_New.mat'], 'n');

%% 2.������������
fij = cell(e_frame-1,1);
fid = cell(e_frame-1,1);
fit = cell(e_frame-1,1);
fiv = cell(e_frame-1,1);
fsj = cell(e_frame,1);
fmj = cell(e_frame,1);


for t = s_frame:e_frame-1
    %  t�е�m��ǰ���еĵ�width��ϸ��   
    fij{t} = binvar(n(t), 4, 'full'); %%fij��������
    fit{t} = binvar(n(t), 1, 'full');   %%��ʧ
    fid{t} = binvar(n(t), 6, 'full');   %%ĸϸ��
    fiv{t} = binvar(n(t), 6, 'full');   %%����
end

for t = s_frame+1:e_frame
    fmj{t} = binvar(n(t), 6,  'full');   %%�ں� ##ע������������t * t-1�����������ǣ�t * t+1��
    fsj{t} = binvar(n(t), 1, 'full');   %%���� 
end

% Ŀ�꺯�� object_f
% w = ones(42,1);
% ============================================
% �¸Ķ���ʹ���������ֻ���� fai_x_z �� sum_cost����ϵĹ��̷ŵ���������
% object_f = w'* fai_x_z + sum_cost;
% toc;

end

function sum_f = cell_dot_mutil( feature, z )
%
% f Ϊ��������cell	n(t)*n(t+1)
% g Ϊ���̱�������
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
sum_f = 0;
ss = numel(feature);
% feature1 = cellfun(@(x)x', feature, 'un',0);
feature2 = reshape(feature, ss, 1);
z1 = reshape(z, ss, 1);
f_z = arrayfun(@(x)feature2{x}*z1(x), 1:ss, 'un',0);
for x=1:ss
    sum_f = sum_f + f_z{x};
end

end











