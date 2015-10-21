%% 函数3：去除长度小于l的edgelist
function edgelist = Remove_Small_Edgelist( edgelist, l)

k=1;
while(k<=numel(edgelist))
    % 正常情况下edgelist首尾相接，如不相接，则手动接上
    if edgelist{k}(1,:)~= edgelist{k}(end,:)
        edgelist{k} = [edgelist{k};edgelist{k}(1,:)];
    end
    while size(edgelist{k},1) < l     % 去除边际小于20的前景
        edgelist(k)=[];
        if k>numel(edgelist)
            break;
        end
    end
    k=k+1;
end


end