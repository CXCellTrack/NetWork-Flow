%% 函数6：压制凹点周边的像素（原先在外部，现在移入）
function [ out ] = CX_sup_aodian(num_ed, ind, sup_distance)  % 输入凹点在edgelist中的编号

% ind=[2 5 8 11 12 340 345 358 359];
% ind即为 ind_In_edgelist    sup_distance为压制长度   num_ed用来做循环压制
%
if numel(ind)==1
    out=ind;
else
    j=1;
    while(j<numel(ind))
        while(ind(j+1)-ind(j)+1<sup_distance)
            ind(j+1)=[];
            if j==numel(ind)
                break;
            end
        end 
        j=j+1;      
    end
    tmp=numel(ind);
    while(numel(ind)==tmp)
        while(ind(1)+num_ed-ind(j-1)+1<sup_distance)           %%all是整个圆周上的点数目
            ind(j-1)=[];
            j=j-1;
        end
        tmp=tmp-1;
    end
    out=ind;
end

end