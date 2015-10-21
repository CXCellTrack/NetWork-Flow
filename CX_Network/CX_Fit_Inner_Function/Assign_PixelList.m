%% ����4������������С���϶Ρ���Ӧ��pixellist
function [ r flag_zuhe ] = Assign_PixelList( r, i, na )

% ���������������  
cp_x=1:na;                  % e.g. 3�����㣬[1,2,3]
cp_sum = zeros(na+1,1);
for jj=1:na+1
    if jj==1
        cp_sum(jj) = 0;
    else
        % e.g. cp_sum ����Ϊ [ 0, C31, C32+C31, C33+C32+C31 ]
        cp_sum(jj) = cp_sum(jj-1) + nchoosek(na,jj-1);
    end
end

if na>=5 % ������Ŀ����ʱ�����������ڸ��ӣ����ֻѡ�����22��ϵ����
    cp_sum = cp_sum(1:3);
    cp_x = 1:2;
end
% e.g. flag_zuhe = false(7,3)
flag_zuhe = false(cp_sum(end), na);
for jj=cp_x
    cp=combntns(1:na,jj);
    for kk=1:size(cp,1)
        flag_zuhe( cp_sum(jj) + kk, cp(kk,:) )=1;
    end
end
% ���������e.g. flag_zuhe = [ 1 0 0;0 1 0;0 0 1 ... ;1 1 1]
% �������˶�������ϣ�1�����λ���ϵ��Ƕ��������������
% ���ɸ������������
n_zuhe = cp_sum(end) - cp_sum(2);
rtemp = cell(1, n_zuhe);
for h=1:n_zuhe
    zuhe = find(flag_zuhe(h+cp_sum(2),:));
    for ii=1:numel(zuhe)
        rtemp{h} = [ rtemp{h}(1:end-1,:); r{i,zuhe(ii)} ];
    end
end
        
% �Ѷ����ϵĽ����r������ϣ�����1,2,3�������12,13,23,123
r(i,(na+1):(na+n_zuhe)) = rtemp;


end