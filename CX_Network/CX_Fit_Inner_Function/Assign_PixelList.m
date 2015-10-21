%% 函数4：计算多段拟合中“组合段”对应的pixellist
function [ r flag_zuhe ] = Assign_PixelList( r, i, na )

% 生成排列组合向量  
cp_x=1:na;                  % e.g. 3个凹点，[1,2,3]
cp_sum = zeros(na+1,1);
for jj=1:na+1
    if jj==1
        cp_sum(jj) = 0;
    else
        % e.g. cp_sum 最终为 [ 0, C31, C32+C31, C33+C32+C31 ]
        cp_sum(jj) = cp_sum(jj-1) + nchoosek(na,jj-1);
    end
end

if na>=5 % 凹点数目过多时，组合情况过于复杂，因此只选择最多22组合的情况
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
% 上面做完后，e.g. flag_zuhe = [ 1 0 0;0 1 0;0 0 1 ... ;1 1 1]
% 即生成了二进制组合，1代表该位置上的那段轮廓参与了拟合
% 生成各段轮廓的组合
n_zuhe = cp_sum(end) - cp_sum(2);
rtemp = cell(1, n_zuhe);
for h=1:n_zuhe
    zuhe = find(flag_zuhe(h+cp_sum(2),:));
    for ii=1:numel(zuhe)
        rtemp{h} = [ rtemp{h}(1:end-1,:); r{i,zuhe(ii)} ];
    end
end
        
% 把多段组合的结果在r后面接上，比如1,2,3，则接上12,13,23,123
r(i,(na+1):(na+n_zuhe)) = rtemp;


end