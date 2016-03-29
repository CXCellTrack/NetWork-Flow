function TRAmea(flowvars_path, dataset)

[ ~, trackpath ] = getpath( dataset );

disp(datestr(now,31));
disp('生成标准格式以供评估...');

% 使用这个近似的投射效果还不错（针对新软件漏检(FN)惩罚10，虚景(FP)惩罚1的特点）
CXSL_Change_to_AOGMM_without_merge( dataset, flowvars_path );

% system(['"%EvaluationSoftware%\SEGMeasure.exe" ', trackpath(1:end-14), ' ', trackpath(end-12:end-11)]);
system(['"%EvaluationSoftware%\TRAMeasure.exe" ', trackpath(1:end-14), ' ', trackpath(end-12:end-11)]);
system(['C:\Users\Administrator\Desktop\AOGMMeasure\AOGMMeasure.exe ',... 
    trackpath(1:end-11),'_GT\TRA ', trackpath(1:end-11),'_RES', ' 5 10 1 1 1.5 1']); % 这个数越小越好？
      
% CXSL_Change_to_AOGMM( flowvars_path ); %
% 这个可以精确分离merge中的2个，将merge和split都变为move，但情况复杂时容易报错
%（比如2个merge后再和别的merge，或是错误的merge导致一直以为是2个