%% ����3��ȥ������С��l��edgelist
function edgelist = Remove_Small_Edgelist( edgelist, l)

k=1;
while(k<=numel(edgelist))
    % ���������edgelist��β��ӣ��粻��ӣ����ֶ�����
    if edgelist{k}(1,:)~= edgelist{k}(end,:)
        edgelist{k} = [edgelist{k};edgelist{k}(1,:)];
    end
    while size(edgelist{k},1) < l     % ȥ���߼�С��20��ǰ��
        edgelist(k)=[];
        if k>numel(edgelist)
            break;
        end
    end
    k=k+1;
end


end