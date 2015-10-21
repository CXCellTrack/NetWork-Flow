function [ segpath trackpath ] = getpath( dataset )

% 返回要寻找的环境变量 pathname 的值
if strcmp(dataset,'training')
    segpath = 'trainsegpath';
    trackpath = 'traintrackpath';
elseif strcmp(dataset,'competition')
    segpath = 'testsegpath';
    trackpath = 'testtrackpath';  
else
    error('输入变量必须为training或competition！');
end

[status1, result1] = system(['set ', segpath]); 
[status2, result2] = system(['set ', trackpath]); 

if ~status1 && ~status2
    segpath = mydeblank( result1((2+numel(segpath)):end) );
    trackpath = mydeblank( result2((2+numel(trackpath)):end) );
else
    error(['没有找到环境变量',segpath,'和',trackpath]);
end