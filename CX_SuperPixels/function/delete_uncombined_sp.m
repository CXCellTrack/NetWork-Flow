function [ label_zuhe flag_zuhe ] = delete_uncombined_sp( fore, label_nearby, label_zuhe, flag_zuhe, bsp_stats, new_labels )

% sp 代表 superpixels
% fore：前景中的所有basic sp
% label_zuhe：全组合，共2^n-1种
% label_nearby：每个label相邻的label集合

% 最终删除不相邻的组合方式
% ======================================= %


%% A方法：使用图来判断组合是否成立（速度稍慢（0.059秒），但保证正确）
% tic
for i_csp=1:numel(label_zuhe) % 对fore中的所有superpixel迭代
    % 当前sp
    csp = label_zuhe{i_csp}; 
    % 绘图底板
    bb = zeros(size(new_labels)); 
    for bsp=csp
        pl = bsp_stats(bsp).PixelList;
        for i_xy=1:size(pl,1) % 将bsp画在图上，再检查图的连通性
            bb(pl(i_xy,1), pl(i_xy,2)) = 1;
        end
    end
    [~, nlabel] = bwlabel(bb, 4);
    
    if numel(csp)==1 % 对于bsp，断言其只有一块，否则就是前面出错了
        assert(nlabel==1)
    end
    if nlabel~=1
        label_zuhe{i_csp} = [];
    end
end

ind_save = ~isemptycell(label_zuhe);
label_zuhe = label_zuhe(ind_save);
flag_zuhe = flag_zuhe(ind_save,:);
% toc

if 0
%% B方法：不需要画图，使用相邻关系来判断组合是否成立（速度快，但暂时有问题）% 2015.12.10 修正了这一方法（0.033秒）
tic
% 组合中有些是不相邻的，要去掉
for bsp=fore % 对fore中的basic_superpixel迭代
    % bsp与它的邻域构成的集合
    S_nb = [label_nearby{bsp}, bsp]; 
    for i_csp=numel(fore)+1:numel(label_zuhe) % 对所有的combine_superpixel迭代(basic_sp不需要检查连通性)
        csp = label_zuhe{i_csp};
        if any(csp==bsp) % 如果组合中包含bsp
            neighbor_in_csp = intersect(S_nb, csp); % csp中包含的bsp和它的邻居
            rest = mysetdiff(csp, S_nb); % 剩下的部分
            if isempty(rest) % 如果没剩下，说明可行，直接跳过
                continue;
            end
            % 如果S_nb中所有的bsp都和rest不相连，则说明rest为单独一块，组合失败
            delete_flag = true;
            for kk=neighbor_in_csp
                tmp = arrayfun(@(x) any(rest==x),label_nearby{kk},'un',1);
                if any(tmp==1) % 只要bsp及其邻居有一个能连接到rest，则说明联通
                    delete_flag = false;
                    break;
                end
            end
            if delete_flag
                label_zuhe{i_csp} = [];
            end
        end
    end
end

ind_save = ~isemptycell(label_zuhe);
label_zuhe = label_zuhe(ind_save);
flag_zuhe = flag_zuhe(ind_save,:);
toc

end
