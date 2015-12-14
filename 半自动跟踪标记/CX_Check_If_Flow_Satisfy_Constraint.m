% ======================================================================= %
% 2015.8.27
% �ڵõ�һ�����̱����Ժ�Fij Fid Fit Fiv Fsj Fmj
% 1��������̱�������cplex���õ��ģ�������������������Լ��
% 2��������̱��������˹�����پ�����ش���õ��ģ���һ���������Լ������Ҫ����ڳ����غ�Լ����
%
% ���ű������þ��Ǽ��õ������̱����Ƿ������Ҫ��Լ������
% ������Ļ������ CX_Visualize_Track_Pair_New ������ɫʱ�ᱨ�������Ҫ�ڴ˴����
% ======================================================================= %

clear;close all;
% ����GT���̱�������Ҫ��������  ����pre_data
if 0
    dataset = 'competition';
else
    dataset = 'training';
end
[ ~, trackpath ] = getpath( dataset );

% ʹ��ȫ�ֱ���
global Fij Fit Fid Fiv Fmj Fsj;
global conflict_fij conflict_pair_last_xy conflict_pair_next_xy n;

load([ trackpath, '\GT\GT_after_hand_tune\GT_Flow_Variables_New.mat']);
load([ trackpath, '\Pair\Pre_data_New.mat'], 'conflict_fij','conflict_pair_last_xy','conflict_pair_next_xy','n','Ellipse','SuperPixel');
frame = numel(Fsj);
frame = 65;

%% ��ʼ��飺��Ҫ����ڳ����غ�

for t=2:frame-1
    disp(['  ����',num2str(t),'֡...']);
    for j=1:n(t)
        % ���ú������j����������
        [ eventIn eventOut ] = CX_CheckInOut( t, j );
        
        sum_enter = sum(eventIn);
        sum_leave = sum(eventOut);
        
        % ��ÿ����Բ����ڳ�����������ж�
        if sum_enter==sum_leave && sum_enter<=1
            % ��ں�=���ںͣ��Ҷ�Ϊ0��Ϊ1����������
            continue;
        else
            % ������������磬��˫��ڵ����ڵ�
            disp(eventIn);
            disp(eventOut);
            
            if exist('SuperPixel','var')
                error(['��', num2str(t), '֡���Ϊ', num2str(j), '�ĳ��������Ϊ',num2str(sum_enter),'��������Ϊ',num2str(sum_leave),...
                ';  ���꣺',num2str(SuperPixel{t}{j}.centroid),';  BSP: ',num2str(SuperPixel{t}{j}.label)]);
            else
                error(['��', num2str(t), '֡���Ϊ', num2str(j), '����Բ���Ϊ',num2str(sum_enter),'��������Ϊ',num2str(sum_leave),...
                    '��  ���꣺',num2str([Ellipse{t}{j}.x0, Ellipse{t}{j}.y0])]);
            end
        end
        
    end
end
disp('  ��Ҫ�������̱�����������غ�Լ��');










