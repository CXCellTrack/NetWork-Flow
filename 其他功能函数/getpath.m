function [ segpath trackpath ] = getpath( dataset )

% ����ҪѰ�ҵĻ������� pathname ��ֵ
if strcmp(dataset,'training')
    segpath = 'trainsegpath';
    trackpath = 'traintrackpath';
elseif strcmp(dataset,'competition')
    segpath = 'testsegpath';
    trackpath = 'testtrackpath';  
else
    error('�����������Ϊtraining��competition��');
end

[status1, result1] = system(['set ', segpath]); 
[status2, result2] = system(['set ', trackpath]); 

if ~status1 && ~status2
    segpath = mydeblank( result1((2+numel(segpath)):end) );
    trackpath = mydeblank( result2((2+numel(trackpath)):end) );
else
    error(['û���ҵ���������',segpath,'��',trackpath]);
end