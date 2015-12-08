function add_into_cluster( label_nearby, row)
%
% 2015.12.8
% 递归函数，用于将所有的邻居都放到一个前景中
%
global cluster

% 单独前景直接添加即可
if isempty(label_nearby{row})
    cluster = [cluster, row];
    return
end
    
for i=label_nearby{row}
    % 先把本行添加进去
    if ~any(cluster==row)
        cluster = [cluster, row];
    end
    % 如i不在c内，则添加，并对i行进行递归调用
    if ~any(cluster==i)
        cluster = [cluster, i];
        add_into_cluster( label_nearby, i);
    end 
    
end