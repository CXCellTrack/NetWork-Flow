function [ seg, ed, flag_aodian, hasError ] = Analyse_aodian( seg, ed, max_aodian )

hasError = 0;

n_seg = size(seg, 1);
% 绘制分段图（只是为了debug方便，并无实际用处）
for j=1:n_seg-1                                    % x向左增大，y向下增大          
    x_j = seg(j,2); y_j = seg(j,1);
    x_k = seg(j+1,2); y_k=seg(j+1,1);
    % 连线，这部分可以不用，后续replot就不绘制了
    plot([x_j x_k],[y_j y_k], 'g', 'LineWidth', 1.5); 
end

%% 1.判断线性分段的每个端点凹点还是凸点（先假设seg为逆时针排列）

%     1 <--- end-1
%    / 
%   < 
%  2 ---> 3 ...
%
% 上图为示意图，∠123为凸角时 det([v1;v2])>0

flag_aodian = false(1, n_seg);          % 预定义判断凹点存在的标识符
for j=1:n_seg-1    
    if j==1 % 第一个点需要特殊处理
        v1=seg(j,:)-seg(n_seg-1,:);       % 前一点到当前点向量
        v2=seg(j+1,:)-seg(j,:);      % 当前点到后一点向量       
    else % 其他点
        v1=seg(j,:)-seg(j-1,:);     
        v2=seg(j+1,:)-seg(j,:);     % edgelist 是逆时针排列么？？
    end
    % 叉乘后第三个向量的大小和方向信息
    v_multi = det([v1;v2]);                 

    if v_multi<0    % 说明为凹点
        flag_aodian(j)= 1;    % 凹点记录
    end              
end

%% 2.判断凹点数目，过多则说明有错误（只能用这种办法来修复 edgelink 函数的bug）

reach_max_aodian = 1; % 确定是否存在连续凹点的标示符
for j=1:n_seg-(max_aodian-1)
    % 连续出现4个凹点，说明 edgelist 并非按正常逆时针在走了
    for ia=0:max_aodian-1
        if flag_aodian(j+ia)==0 % 连续几个点中，只要有一个非凹，就认为不存在连续凹点的错误情况
            reach_max_aodian = 0;
        end
    end
    
    if reach_max_aodian
        seg = fliplr(seg')';  % 将序列反转，变为顺时针
        ed = fliplr(ed')';
        flag_aodian = [ ~flag_aodian(1), ~fliplr(flag_aodian(2:end)) ];
        break;
    end
end

for j=1:n_seg-4
    % 反转后还是连续出现4个凹点，说明出现了8字形
    if flag_aodian(j) && flag_aodian(j+1) && flag_aodian(j+2) && flag_aodian(j+3) && flag_aodian(j+4)
        hasError = 1;
        return;
    end
end
        
% 如果上述2步检查都无误，则说明无错误
% ------------------------- 可删 ------------------------------- %
for j=1:numel(flag_aodian)
    if flag_aodian(j)
        plot(seg(j,2), seg(j,1), 'bs', 'LineWidth', 1.5); % 凹点画一个X做标记（可以不画）
    end
end
% --------------------------------------------------------------- %


% %% 1.检查 edgelist 是否为逆时针排列，如果不是，将其转化为逆时针
% 
% %     1 ---> end-1
% %    / 
% %   < 
% %  2 ---> 3 ...
% %
% % 上图为示意图，逆时针排列时 det([p21;pn1])>0
% p21 = seg(2,:)-seg(1,:);
% pn1 = seg(end-1,:)-seg(1,:);
% % det代表叉乘运算
% if det([p21;pn1])<0     % 小于0说明seg是顺时针的，需要翻一下,同时edgelist也需要反转
%     seg=fliplr(seg')';  % (2015.4.29)修正了前景相连时edgelink函数产生的bug
%     ed=fliplr(ed')';
% end
% % 上述修正后将所有 edgelist 和 seglist 都变成了逆时针











