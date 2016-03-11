function  e = CX_Ellipse_Optimal( e )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% 下面一行可用于debug
% e = Ellipse;
frame = numel(e);

for t=1:frame 
    disp(['优化筛选第',num2str(t),'帧...']);
    flag = ~isemptycell(e{t});
    height = size(flag,1);
    for hh=1:height        %%每个前景
        
        width = sum(flag(hh,:));  %%计算每行的cell数目，即前景的拟合数目
        
        %% 第一步：删去status为0的假说
        %########### 删去status为0的假说 ##########
		tmp_e = e{t}(hh,:);
		for i=1:width
			if tmp_e{i}.status==0
				tmp_e{i} = [];
			end
		end
		sts_flag = ~isemptycell(tmp_e);
        new_width = sum(sts_flag); % 删除后的长度
        if new_width==0 % 整行都删没了
            e{t}(hh,1) = {[]};
            continue
        end  
        e{t}(hh,1:new_width) = tmp_e(sts_flag);
        e{t}(hh,new_width+1:end) = {[]};   

        %% 第二步：在多凹点情况下椭圆假说太多，删除一些拟合程度较差的
         %######## 打分函数，很重要！！######## 
        tmp_e = tmp_e(sts_flag);
        % B. 直接在cell基础上做排序,方便
        tmp_fitness = zeros(1,new_width);
        for i=1:new_width
            % 打分函数暂定为 hd/sqrt(num_pixels)
            tmp_fitness(i) = tmp_e{i}.hd/sqrt(tmp_e{i}.num_pixels);
        end
        [~, ind_fit] = sort(tmp_fitness);
        sort_e = tmp_e;
        for i=1:new_width
            sort_e(i) = tmp_e(ind_fit(i));
        end
        % 排序完的一个前景可以进行删减，再返回原cell
        if new_width >15
            num_remain = 15; % 假说大于10个则保留拟合度最好的10个，否则全部保留
        else
            num_remain = new_width; % 定义保留的个数
        end
        % 给e赋值
        e{t}(hh, 1:num_remain) = sort_e(1:num_remain);
        e{t}(hh, num_remain+1:end) = {[]};
        
        %% 第三步：寻找相似椭圆 （2015.6.27修改了一个bug，将下面的 new_width 改正为 num_remain ）
        for ind_1=1:num_remain-1
            for ind_2=ind_1+1:num_remain
                e1=e{t}{hh,ind_1};
                e2=e{t}{hh,ind_2};
                % 如果椭圆无交集，直接进入下一轮循环
                if max([abs(e1.x0-e2.x0), abs(e1.y0-e2.y0)]) > e1.a + e2.a
                    continue;
                end

%                 e1_4_c = struct2array(e1);
%                 e2_4_c = struct2array(e2);   
                % 新发现这个函数有点问题，无法正确计算出来，先不用了
%                 overlap = CX_calculate_overlap_C(e1_4_c,e2_4_c);%%使用C编写的重叠度矩阵，72秒，改为计算面积，变成17秒

                overlap = CX_calculate_overlap(e1,e2);   % CX_calculate_overlap 是matlab版的计算过程，166秒
                if overlap >= 0.75    %%重叠度定义为0.9
%                     belta_i = e1.num_pixels;
%                     belta_j = e2.num_pixels;   %%直接比较拟合边的长度
                    e{t}{hh,ind_2}.similar_delete = 1; % 因为e2排在后面，说明拟合程度较差，因此将其删除
                    
                if 0 % 注释掉的方法
                  %####  判断圆心角的方法速度太慢
    %             %%判定为同一个时则比较圆心角大小
    %                 rr=r{n,j}*[0 1;1 0];   %%rr(y,x)需要换列才能变成rr(x,y)
    %                 v1=rr(1,:)-[e2.x0,e2.y0];
    %                 v2=rr(round(end/2),:)-[e2.x0,e2.y0];
    %                 v3=rr(end,:)-[e2.x0,e2.y0];
    %                 theta1=acosd(v1*v2'/sqrt(v1*v1'*v2*v2'));
    %                 theta2=acosd(v2*v3'/sqrt(v2*v2'*v3*v3'));
    %                 belta_j=theta1+theta2;
    %                 temp_e{n,i}.theta=belta_i; 
    
%                     if belta_i>=belta_j
%                         e{t}{n,j}.similar_delete=1;   %%将拟合边长度小的那个similar_delete=1，在后续就不计算了
%                     else
%                         e{t}{n,i}.similar_delete=1;
%                     end
                end
                
                end
            end
        end
        
        %% 第四步：删除相似椭圆，更新 num_hypoth 用来判断该前景是否为单假说前景
        % --------------------------------------------------------------- %
        % 下面这部分很重要，不更新 num_hypoth 的话会导致后面人工标记时出错
        % 并且需要利用 num_hypoth 判断该椭圆是否位于多假说前景中
        % --------------------------------------------------------------- %
        % 添加删去 similar_delete=1 的椭圆，并重新统计 num_hypoth  2015.6.28
		tmp_e = e{t}(hh,:);
		for i=1:num_remain
			if tmp_e{i}.similar_delete==1
				tmp_e{i} = [];
			end
		end
		simd_flag = ~isemptycell(tmp_e);
        num_remain = sum(simd_flag); % 删除后的长度
		e{t}(hh,1:num_remain) = tmp_e(simd_flag);
        e{t}(hh,num_remain+1:end) = {[]};

        % 需更新前景内所有椭圆的 num_hypoth
        for ii=1:num_remain
            e{t}{hh,ii}.num_hypoth = num_remain;
        end
        
        %% 第五步：按照 sum（flag_combine）重新排序，将单段拟合的放在前面
        % --------------------------------------------------------------- %
        % 总结 2015.7.7 ---
        
        % 运行到这里时，每个前景中的 status 为0和 similar_delete 为1的椭圆已经
        % 删去，并且按照拟合程度取前10位进行排序（大于10个的截取10个，小于10个的取原数目）
        %
        % 然而在后面生成矛盾假说时，因为是从前往后取 numel(flag_combine) 个单1椭圆作为基组
        % 例如三段中的 100 010 001，并扩展出3个矛盾集合
        % 
        % 因此： 在下面的代码中需要对每一行的椭圆进行重新排序，使得单点在前面。
        % *******（注意，在这个方法没实现时使用的是 22矛盾的假说，容易造成约束数量过多，影响计算速度） ******
        
        tmp_e = e{t}(hh, 1:num_remain);
        sum_flag_combine = zeros(1, num_remain);
        % 求出 flag_combine 的和
        for ind_e=1:num_remain
            sum_flag_combine(ind_e) = sum( tmp_e{ind_e}.flag_combine );
        end
        % 将和小的放在前面，即单段放前面
        [ ~, new_ind ] = sort(sum_flag_combine);
        % -------- %
        % 将值赋给e，过程结束
        e{t}(hh, 1:num_remain) = tmp_e(new_ind);
        e{t}(hh, num_remain+1:end) = {[]};
        % --------------------------------------------------------------- %
        
    end
    %% 其他优化属性（暂无）
      
end

end % end function


%% 已经改写为C版本 下面这个matlab版本不用了（后来发现c版本有问题，因此仍使用下面）

function [ overlap ] = CX_calculate_overlap(e1,e2)

        c1=sqrt(e1.a^2-e1.b^2);   %%焦距
        cc1=cosd(e1.alpha);
        ss1=sind(e1.alpha);
        c11=[c1*cc1+e1.x0,c1*ss1+e1.y0];
        c12=[-c1*cc1+e1.x0,-c1*ss1+e1.y0];
        
        count1 = pi*e1.a*e1.b;
        count2 = pi*e2.a*e2.b;
        count=0;
        
        c2=sqrt(e2.a^2-e2.b^2);   %%焦距
        cc2=cosd(e2.alpha);
        ss2=sind(e2.alpha);
        c21=[c2*cc2+e2.x0,c2*ss2+e2.y0];
        c22=[-c2*cc2+e2.x0,-c2*ss2+e2.y0];

        for x=e2.x0-e2.a : e2.x0+e2.a
            for y=e2.y0-e2.a : e2.y0+e2.a
                if sqrt((x-c21(1))^2+(y-c21(2))^2)+sqrt((x-c22(1))^2+(y-c22(2))^2)<=2*e2.a
                    if sqrt((x-c11(1))^2+(y-c11(2))^2)+sqrt((x-c12(1))^2+(y-c12(2))^2)<=2*e1.a
                        count=count+1;
                    end
                end
            end
        end
        overlap=count/(count1+count2-count);
end

