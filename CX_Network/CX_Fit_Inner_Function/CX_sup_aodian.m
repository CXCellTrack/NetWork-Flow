%% ����6��ѹ�ư����ܱߵ����أ�ԭ�����ⲿ���������룩
function [ out ] = CX_sup_aodian(num_ed, ind, sup_distance)  % ���밼����edgelist�еı��

% ind=[2 5 8 11 12 340 345 358 359];
% ind��Ϊ ind_In_edgelist    sup_distanceΪѹ�Ƴ���   num_ed������ѭ��ѹ��
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
        while(ind(1)+num_ed-ind(j-1)+1<sup_distance)           %%all������Բ���ϵĵ���Ŀ
            ind(j-1)=[];
            j=j-1;
        end
        tmp=tmp-1;
    end
    out=ind;
end

end