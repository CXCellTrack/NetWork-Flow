function [ label_zuhe flag_zuhe ] = delete_combined_sp( fore, label_nearby, label_zuhe, flag_zuhe )

% sp 代表 superpixels
% fore：前景中的所有basic sp
% label_zuhe：全组合，共2^n-1种
% label_nearby：每个label相邻的label集合

% 最终删除不相邻的组合方式
% ======================================= %

% 组合中有些是不相邻的，要去掉
for bsp=fore % 对fore中的basic_superpixel迭代
    for i_csp=numel(fore)+1:numel(label_zuhe) % 对所有的combine_superpixel迭代
        csp = label_zuhe{i_csp};
        if any(csp==bsp) % 如果组合中包含bsp
            rest = mysetdiff(csp,bsp); % 除了bsp剩下的的与bsp的邻域做交集
            % A方法：使用集合运算速度稍微慢一点
            % 如果相交为空，说明rest和bsp都不相邻，则删除这种假说
%                 if isempty( intersect(rest, label_nearby{bsp}) ) 
%                     label_zuhe{i_csp} = [];
%                 end
            % B方法：arraufun逐个比较，速度快
            tmp = arrayfun(@(x) any(rest==x),label_nearby{bsp},'un',1);
            if all(tmp==0)
                label_zuhe{i_csp} = [];
            end
        end
    end
end

ind_save = ~isemptycell(label_zuhe);
label_zuhe = label_zuhe(ind_save);
flag_zuhe = flag_zuhe(ind_save,:);

