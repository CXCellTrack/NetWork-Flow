function add_into_cluster( label_nearby, row)
%
% 2015.12.8
% �ݹ麯�������ڽ����е��ھӶ��ŵ�һ��ǰ����
%
global cluster

% ����ǰ��ֱ����Ӽ���
if isempty(label_nearby{row})
    cluster = [cluster, row];
    return
end
    
for i=label_nearby{row}
    % �Ȱѱ�����ӽ�ȥ
    if ~any(cluster==row)
        cluster = [cluster, row];
    end
    % ��i����c�ڣ�����ӣ�����i�н��еݹ����
    if ~any(cluster==i)
        cluster = [cluster, i];
        add_into_cluster( label_nearby, i);
    end 
    
end