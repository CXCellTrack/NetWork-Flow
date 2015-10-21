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
if 1
    dataset = 'competition';
else
    dataset = 'training';
end
[ ~, trackpath ] = getpath( dataset );

load([ trackpath, '\GT\GT_after_hand_tune\GT_Flow_Variables_New.mat']);
load([ trackpath, '\Pair\Pre_data_New.mat']);
frame = numel(Fmj);
% frame = 40;

%% ��ʼ��飺��Ҫ����ڳ����غ�

for t=2:frame-1
    disp(['  ����',num2str(t),'֡...']);
    for j=1:n(t)
        sum_fij = 0;
        sum_fid = 0;
        sum_fiv = 0;
        sum_fmj = 0;

        if t>1 % ��һ֡û�����
            % sum_fij Ϊ����Ǩ�Ƶ� j ��fij֮�ͣ���ںͣ������ںͿ����� sum(fij{t}(j,:)) ��ʾ 
            for ind=1:size(conflict_fij{t-1}{j}, 1)
                sum_fij = sum_fij + Fij{t-1}( conflict_fij{t-1}{j}(ind,1), conflict_fij{t-1}{j}(ind,2) );
            end

            % sum_fid Ϊ���з��ѵ����� j �� pair �� fid ֮�ͣ���ڣ�
            for ind=1:numel(conflict_pair_next_xy{t}{j})/2
                sum_fid = sum_fid + Fid{t-1}( conflict_pair_next_xy{t}{j}(ind,1), conflict_pair_next_xy{t}{j}(ind,2) );
                sum_fiv = sum_fiv + Fiv{t-1}( conflict_pair_next_xy{t}{j}(ind,1), conflict_pair_next_xy{t}{j}(ind,2) );
            end
            
            eventIn(1) = sum_fij;
            eventIn(2) = Fsj{t}(j);
            eventIn(3) = sum(Fmj{t}(j,:));
            eventIn(4) = sum_fid;
            eventIn(5) = sum_fiv;
            sum_enter = sum_fij + Fsj{t}(j) + sum(Fmj{t}(j,:)) + sum_fid + sum_fiv; % j������ں�
        end
        
        if t<frame % ���һ֡û�г���
            % sum_fmj Ϊ���а��� j ���ں� pair �� fmj ֮�ͣ����ڣ�
            for ind=1:numel(conflict_pair_last_xy{t}{j})/2
                sum_fmj = sum_fmj + Fmj{t+1}( conflict_pair_last_xy{t}{j}(ind,1), conflict_pair_last_xy{t}{j}(ind,2) );
            end
            
            eventOut(1) = sum(Fij{t}(j,:));
            eventOut(2) = Fit{t}(j);
            eventOut(3) = sum(Fid{t}(j,:));
            eventOut(4) = sum(Fiv{t}(j,:));
            eventOut(5) = sum_fmj;

            sum_leave = sum(Fij{t}(j,:)) + Fit{t}(j) + sum(Fid{t}(j,:)) + sum(Fiv{t}(j,:)) + sum_fmj; % j���ܳ��ں�
        end
        

        % ��ÿ����Բ����ڳ�����������ж�
        if sum_enter==sum_leave && sum_enter<=1
            % ��ں�=���ںͣ��Ҷ�Ϊ0��Ϊ1����������
            continue;
        else
            % ������������磬��˫��ڵ����ڵ�
            disp(eventIn);
            disp(eventOut);
            error(['��', num2str(t), '֡���Ϊ', num2str(j), '����Բ���Ϊ',num2str(sum_enter),'��������Ϊ',num2str(sum_leave),...
                '��  ���꣺',num2str([Ellipse{t}{j}.x0, Ellipse{t}{j}.y0])]);
        end
        
    end
end
disp('  ��Ҫ�������̱�����������غ�Լ��');










