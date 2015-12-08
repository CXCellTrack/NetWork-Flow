function [ tmpe, n ] = CX_Cut_Ellipse( Ellipse )

% ����������ڳ�ȥraw_ellipse�еĿ�cell������ֱcell

tmpe = cell(size(Ellipse)); % ��ʱ�����Բcell
n = zeros(size(Ellipse)); % ÿ֡�ڵ���Բ����

for i=1:numel(Ellipse)
    width = size(Ellipse{i},2);
    flag = ~isemptycell(Ellipse{i}); % �������ⲿ����isemptycell
    
    for j=1:size(flag,1)
        for k=1:sum(flag(j,:))
            n(i) = n(i) + 1;
            tmpe{i}{n(i),1} = Ellipse{i}{j,k};
        end
    end
end