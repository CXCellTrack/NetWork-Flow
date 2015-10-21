function e_used = find_ellipse_used_in_GT()


dataset = 'competition';
[ ~, trackpath ] = getpath( dataset );
% ���� Ellipse ������
load([ trackpath, '\Pair\Pre_data_New.mat']);
% �������յ����̱���
track_data_addr = [ trackpath, '\GT\GT_after_hand_tune\GT_Flow_Variables_New.mat'];
load(track_data_addr);

%% ��ʼ���ÿ����Բ��˵�Ƿ����������
e_used  = [];
frame = numel(Fmj);

for t=1:frame
    ind_used = 0;
    for j=1:n(t)

        sum_fij = 0;
        sum_fid = 0;
        sum_fiv = 0;
        sum_fmj = 0;
        sum_enter = 0;
        sum_leave = 0;

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
            
            sum_enter = sum_fij + Fsj{t}(j) + sum(Fmj{t}(j,:)) + sum_fid + sum_fiv; % j������ں�
        end
        
        if t<frame % ���һ֡û�г���
            % sum_fmj Ϊ���а��� j ���ں� pair �� fmj ֮�ͣ����ڣ�
            for ind=1:numel(conflict_pair_last_xy{t}{j})/2
                sum_fmj = sum_fmj + Fmj{t+1}( conflict_pair_last_xy{t}{j}(ind,1), conflict_pair_last_xy{t}{j}(ind,2) );
            end
            
            sum_leave = sum(Fij{t}(j,:)) + Fit{t}(j) + sum(Fid{t}(j,:)) + sum(Fiv{t}(j,:)) + sum_fmj; % j���ܳ��ں�
        end

        
        % ��ÿ����Բ����ڳ�����������ж�
        if sum_enter || sum_leave
            ind_used = ind_used + 1;
            e_used(t,ind_used) = j; % �����¼��ÿ֡����Բʹ�����
        end
        
    end
end










