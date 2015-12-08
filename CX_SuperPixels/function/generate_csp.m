function [ label_zuhe, flag_zuhe ] = generate_csp( fore )

% 生成所有的组合情况 combined superpixel
na = numel(fore);
cp_x=1:na;                  % e.g. 3个超像素，[1,2,3]
cp_sum = zeros(na+1,1);
for jj=1:na+1
    if jj==1
        cp_sum(jj) = 0;
    else
        % e.g. cp_sum 最终为 [ 0, C31, C32+C31, C33+C32+C31 ]
        cp_sum(jj) = cp_sum(jj-1) + nchoosek(na,jj-1);
    end
end

if na>5 % 超像素数目过多时，组合情况过于复杂，因此只选择最多33组合的情况
    cp_sum = cp_sum(1:4); % 实际上10个超像素33组合也有175种情况，还是非常多
    cp_x = 1:3;
end

% e.g. flag_zuhe = false(7,3)
flag_zuhe = false(cp_sum(end), na);
label_zuhe = cell(cp_sum(end),1);

for jj=cp_x
    cp=combntns(1:na,jj); % 表示二进制组合的cp
    for kk=1:size(cp,1)
        flag_zuhe( cp_sum(jj) + kk, cp(kk,:) ) = 1; % flag_zuhe用来编写矛盾假说比较方便
    end
    
    cp_label = combntns(fore, jj);
    for kk=1:size(cp_label,1)
        label_zuhe{cp_sum(jj)+kk} = cp_label(kk,:); % label_zuhe更为直观
    end
end