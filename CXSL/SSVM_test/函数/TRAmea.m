function TRAmea(flowvars_path, dataset)

[ ~, trackpath ] = getpath( dataset );

disp(datestr(now,31));
disp('���ɱ�׼��ʽ�Թ�����...');

% ʹ��������Ƶ�Ͷ��Ч����������������©��(FN)�ͷ�10���龰(FP)�ͷ�1���ص㣩
CXSL_Change_to_AOGMM_without_merge( dataset, flowvars_path );

% system(['"%EvaluationSoftware%\SEGMeasure.exe" ', trackpath(1:end-14), ' ', trackpath(end-12:end-11)]);
system(['"%EvaluationSoftware%\TRAMeasure.exe" ', trackpath(1:end-14), ' ', trackpath(end-12:end-11)]);
system(['C:\Users\Administrator\Desktop\AOGMMeasure\AOGMMeasure.exe ',... 
    trackpath(1:end-11),'_GT\TRA ', trackpath(1:end-11),'_RES', ' 5 10 1 1 1.5 1']); % �����ԽСԽ�ã�
      
% CXSL_Change_to_AOGMM( flowvars_path ); %
% ������Ծ�ȷ����merge�е�2������merge��split����Ϊmove�����������ʱ���ױ���
%������2��merge���ٺͱ��merge�����Ǵ����merge����һֱ��Ϊ��2��