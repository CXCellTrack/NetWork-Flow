%% 函数1：划定foi感兴趣的区域
function bw = Remove_Cells_Out_Of_FOI( bw, foi )

if foi==0
    return;
end
% foi = 28;

[ height, width ] = size(bw);
L = bwlabeln(bw,4);
stats = regionprops(L, 'BoundingBox');
for j=1:numel(stats)
    % 算出区域的四个顶点
    tmpbd = stats(j).BoundingBox;
    ltop = tmpbd(1:2);
    lbtm = [ tmpbd(1), tmpbd(2)+tmpbd(4) ];
    rtop = [ tmpbd(1)+tmpbd(3), tmpbd(2) ];
    rbtm = [ tmpbd(1)+tmpbd(3), tmpbd(2)+tmpbd(4) ];
    % 如果四个顶点都在foi外，则将此区域删除
    center = [ width/2, height/2 ]; % 图片中心的坐标
    if min( [abs(ltop(1)-center(1)),abs(lbtm(1)-center(1)),abs(rtop(1)-center(1)),abs(rbtm(1)-center(1))] ) >= center(1)-foi ...
        || min( [abs(ltop(2)-center(2)),abs(lbtm(2)-center(2)),abs(rtop(2)-center(2)),abs(rbtm(2)-center(2))] ) >= center(2)-foi
        bw(L==j)=0;
    end
end

end