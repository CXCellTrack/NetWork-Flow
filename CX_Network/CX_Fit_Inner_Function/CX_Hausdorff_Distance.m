function [ e ] = CX_Hausdorff_Distance( e, r )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%% 凹点数目与椭圆总数关系表
% 凹点数目    椭圆总数
%   1           1
%   2           3
%   3           7
%   4           15
%   5           31
%%

flag = ~isemptycell(e);
% for n=1:size(e,1) % 每个前景
width = sum(flag); % 前景中的假说数目
    
for i=1:width
    if e{i}.status==0 % status 为0 则跳过
        continue;
    end
    e{i} = calculate_distance(e{i} ,r{i}); % 输入当前圆，当前点集，当前全集r{n,width}（似乎用不上）
end

% end
end



function  [ e_n_i ]= calculate_distance( e_n_i ,r_n_i )  %%必须要把值传递出来，不能地址传递

    [ exy ] = ellipse_xy(e_n_i);
    rr = r_n_i*[0 1;1 0];        %%将轮廓矩阵2列交换   r是继承自edgelist,seglist，2个都是xy相反的，需要换位
%     r_all = r_n_width*[0 1;1 0];   %%r_all表示整段轮廓
    [hd ~] = HausdorffDist(rr ,exy, 1);    %%此函数已经更改为计算椭圆包围轮廓的最大最短距离，即贴合度概念
    %% 计算g（cl，e），目前有点问题，先不用
%     complement_rr=setdiff(r_all, rr, 'rows');  %%rr的补集，用于计算g（cl，e）
%     flag= false(numel(complement_rr)/2,1);  %%标志矩阵
%     for j=1:numel(complement_rr)/2
%         xg=complement_rr(j,1);
%         yg=complement_rr(j,2);
%         if sqrt((complement_rr(j)-c11)*(complement_rr(j)-c11)')+sqrt((complement_rr(j)-c12)*(complement_rr(j)-c12)')<=2*ee.a  %%提取rr的补集中位于椭圆内的点
%             flag(j)=1;
%         end
%     end
%     CL=complement_rr(flag,:);
%     if ~isempty(CL)
%         [gd ~] = HausdorffDist(CL,exy,1,'vis');   %%计算g（cl，e）,空矩阵就不算了，加快速度
%     else
%         gd = 0;
%     end
    %%
%     gd=0;
    e_n_i.hd = hd;
    e_n_i.num_pixels = numel(rr)/2; % 2015.5.25
%     e_n_i.shape_yueshu = shape_yueshu;
%     e_n_i.Cl_distance = (hd + gd)*shape_yueshu;  %%距离公式
    
end

function [ exy ]= ellipse_xy(ee)   %%输入椭圆e{i,j}，输出椭圆上的点集和公式后2项,还有2个焦点
    % 返回椭圆周上的点坐标
    c=cosd(ee.alpha);
    s=sind(ee.alpha);
    polar_angle=linspace(0,180,181);
    xq= ee.a*cosd(polar_angle);
    yq= ee.b*sind(polar_angle);
    xn=xq*c-yq*s+ee.x0;
    yn=xq*s+yq*c+ee.y0;
    exy=[xn',yn'];    %%得到椭圆周上的点

%     Ce=pi*sqrt(2*(ee.a^2+ee.b^2));    %%椭圆周长近似公式
%     ec=sqrt(ee.a^2-ee.b^2)/ee.a;   %%偏心   c/a
%     shape_yueshu=Ce*sqrt(1/(1+ec^2));
%     
% 	ccc=sqrt(ee.a^2-ee.b^2);   %%半焦距c
%     c11=[ccc*c+ee.x0 ,ccc*s+ee.y0];
%     c12=[-ccc*c+ee.x0 ,-ccc*s+ee.y0];

    
end
 

    
    
    
