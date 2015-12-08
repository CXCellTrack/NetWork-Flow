function [ label_zuhe, flag_zuhe ] = generate_csp( fore )

% �������е������� combined superpixel
na = numel(fore);
cp_x=1:na;                  % e.g. 3�������أ�[1,2,3]
cp_sum = zeros(na+1,1);
for jj=1:na+1
    if jj==1
        cp_sum(jj) = 0;
    else
        % e.g. cp_sum ����Ϊ [ 0, C31, C32+C31, C33+C32+C31 ]
        cp_sum(jj) = cp_sum(jj-1) + nchoosek(na,jj-1);
    end
end

if na>5 % ��������Ŀ����ʱ�����������ڸ��ӣ����ֻѡ�����33��ϵ����
    cp_sum = cp_sum(1:4); % ʵ����10��������33���Ҳ��175����������Ƿǳ���
    cp_x = 1:3;
end

% e.g. flag_zuhe = false(7,3)
flag_zuhe = false(cp_sum(end), na);
label_zuhe = cell(cp_sum(end),1);

for jj=cp_x
    cp=combntns(1:na,jj); % ��ʾ��������ϵ�cp
    for kk=1:size(cp,1)
        flag_zuhe( cp_sum(jj) + kk, cp(kk,:) ) = 1; % flag_zuhe������дì�ܼ�˵�ȽϷ���
    end
    
    cp_label = combntns(fore, jj);
    for kk=1:size(cp_label,1)
        label_zuhe{cp_sum(jj)+kk} = cp_label(kk,:); % label_zuhe��Ϊֱ��
    end
end