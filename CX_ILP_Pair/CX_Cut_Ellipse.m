function [ tmpe, n ] = CX_Cut_Ellipse( Ellipse )

% 这个函数用于除去raw_ellipse中的空cell，并拉直cell

tmpe = cell(size(Ellipse)); % 临时存放椭圆cell
n = zeros(size(Ellipse)); % 每帧内的椭圆总数

for i=1:numel(Ellipse)
    width = size(Ellipse{i},2);
    flag = ~isemptycell(Ellipse{i}); % 调用了外部函数isemptycell
    
    for j=1:size(flag,1)
        for k=1:sum(flag(j,:))
            n(i) = n(i) + 1;
            tmpe{i}{n(i),1} = Ellipse{i}{j,k};
        end
    end
end