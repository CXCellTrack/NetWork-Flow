% ����ű���������Ŀ¼������ѵ������Կ�ʼǰʹ��һ�μ��ɣ�
clear

%% ѵ��Ŀ¼
dataset = 'training';
[ ~, trackpath ] = getpath( dataset );
path = {};
path{1} = [ trackpath, '\GT'];
path{2} = [ trackpath, '\GT\GT_after_hand_tune'];
path{3} = [ trackpath, '\GT\label_and_e'];
path{4} = [ trackpath, '\Pair'];
path{5} = [ trackpath, '\Pair\���ӻ����ٱ��'];
path{6} = [ trackpath, '\�ṹ��ѧϰ'];
path{7} = [ trackpath, '\�����ͼ'];
path{8} = [ trackpath, '\ѵ�������¼'];
path{9} = [ trackpath, '\ѵ�������¼\BCFWavg_New'];

for i=1:numel(path)
    mkdir(path{i});
end

%% ����Ŀ¼
dataset = 'competition';
[ ~, trackpath ] = getpath( dataset );
path = {};
path{1} = [ trackpath, '\GT'];
path{2} = [ trackpath, '\GT\GT_after_hand_tune'];
% path{3} = [ trackpath, '\GT\label_and_e'];
path{4} = [ trackpath, '\Pair'];
path{5} = [ trackpath, '\Pair\���ӻ����ٱ��'];
path{6} = [ trackpath, '\�ṹ��ѧϰ'];
path{7} = [ trackpath, '\�����ͼ'];
path{8} = [ trackpath, '\���Խ����¼'];
path{9} = [ trackpath, '\���Խ����¼\BCFWavg_New'];


for i=1:numel(path)
    mkdir(path{i});
end