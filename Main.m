% �������ڲ��Լ��ϴӳ�ʼ�ָ�������ո���ͼƬ
% ע�⣬�ű��еĲ���ֻ�ܵ������ļ����޸�

clear;close all;

%% 1��������Բ��˵ raw_ellipse
run('CX_Network');

%% 2���Ż���Բ������ pre_data ������
run('CX_ILP_Pair_Pre_New');

%% 3���������������� feature_plus ������
run('CXSL_Combine_feature_all_New');

%% 4����ѵ���õ���w�����ILP����
run('CXSL_ILP_Using_Best_W_New');

%% 5��������Բ�������ͼ
CX_RePlot_Ellipse( 'competition');

%% 6���������յĸ��ٽ��
run('CX_Visualize_Track_Pair_New');





