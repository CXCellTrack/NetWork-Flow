function [ eventIn, eventOut ] = CX_CheckInOut( t, jj )
%
% 2015.12.9 ����ļ���ôĪ������û�°벿����
% ������

% ʹ��ȫ�ֱ���
global Fij Fit Fid Fiv Fmj Fsj;
global conflict_fij conflict_pair_last_xy conflict_pair_next_xy n;
frame = numel(Fsj);

sum_fij = 0;
sum_fid = 0;
sum_fiv = 0;
sum_fmj = 0;
eventIn = zeros(1,6);
eventOut = zeros(1,6);

if t>1 % ��һ֡û�����
    % sum_fij Ϊ����Ǩ�Ƶ� j ��fij֮�ͣ���ںͣ������ںͿ����� sum(fij{t}(j,:)) ��ʾ 
    for ind=1:size(conflict_fij{t-1}{jj}, 1)
        sum_fij = sum_fij + Fij{t-1}( conflict_fij{t-1}{jj}(ind,1), conflict_fij{t-1}{jj}(ind,2) );
    end

    % sum_fid Ϊ���з��ѵ����� j �� pair �� fid ֮�ͣ���ڣ�
    for ind=1:numel(conflict_pair_next_xy{t}{jj})/2
        sum_fid = sum_fid + Fid{t-1}( conflict_pair_next_xy{t}{jj}(ind,1), conflict_pair_next_xy{t}{jj}(ind,2) );
        sum_fiv = sum_fiv + Fiv{t-1}( conflict_pair_next_xy{t}{jj}(ind,1), conflict_pair_next_xy{t}{jj}(ind,2) );
    end

    eventIn(1) = sum_fij;
    eventIn(2) = 0;
    eventIn(3) = sum_fid;
    eventIn(4) = sum_fiv;
    eventIn(5) = sum(Fmj{t}(jj,:));
    eventIn(6) = Fsj{t}(jj);
    
    sum_enter = sum(eventIn); % j������ں�
end

if t<frame % ���һ֡û�г���
    % sum_fmj Ϊ���а��� j ���ں� pair �� fmj ֮�ͣ����ڣ�
    for ind=1:numel(conflict_pair_last_xy{t}{jj})/2
        sum_fmj = sum_fmj + Fmj{t+1}( conflict_pair_last_xy{t}{jj}(ind,1), conflict_pair_last_xy{t}{jj}(ind,2) );
    end

    eventOut(1) = sum(Fij{t}(jj,:));
    eventOut(2) = Fit{t}(jj);
    eventOut(3) = sum(Fid{t}(jj,:));
    eventOut(4) = sum(Fiv{t}(jj,:));
    eventOut(5) = sum_fmj;
    eventOut(6) = 0;
    sum_leave = sum(eventOut); % j���ܳ��ں�
end











