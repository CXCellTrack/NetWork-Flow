%######################################
%
% 2015.6.5 CX on desk
% 作用：这个函数仅用于被 CX_Network 调用
% 2015.6.5（加入了 FOI 提取，使边缘的判断符合标准 ）
%
%######################################

function [ e ,lunkuo ]= CX_Fit( frame, pic ,iteration_num, tolerance, remove_small )   %%iteration_num是众数滤波次数
    
tic;

    %% 1.读入图片
%     error(nargchk(4, 4, nargin));
    if ischar(pic)
        fprintf('正在处理 %s ...滤波次数 %d\n',pic, iteration_num);  
        im = imread(pic);
    else
        im = pic;
    end
    
    %% 2.二值化
    bw = im2bw(im);
    
    %% 3.设置foi(感兴趣的区域),将部分或全部在区域内的留下，其余的删去 2015.6.5 根据welcome要求制定
%   ? E = 50 for the datasets DIC-C2DH-HeLa, Fluo-C2DL-MSC, Fluo-C3DH-H157, Fluo-N2DHGOWT1,
%   Fluo-N3DH-CE, Fluo-N3DH-CHO, and PhC-C2DH-U373,
%   ? E = 25 for the datasets Fluo-C3DL-MDA231, Fluo-N2DL-HeLa, and PhC-C2DL-PSC,
%   ? E = 0 for the simulated datasets (i.e., Fluo-N2DH-SIM, Fluo-N2DH-SIM+, Fluo-N3DH-SIM, and
%   Fluo-N3DH-SIM+) and for the dataset Fluo-N3DL-DRO.
    foi = 0;
    % 调用内部函数
    bw = Remove_Cells_Out_Of_FOI( bw, foi );

    %% 4.进行形态学上的一些操作
%   % 先清除边际，防止close的时候细胞相连
%   bwbord=imclearborder(bw,4);
%   % 形态学操作容易使原本分开的细胞粘合，因此需谨慎使用   
% 	se=strel('diamond',1);
%   bw_m=imclose(bwbord,se);
    bw_m = bwmorph(bw, 'majority', iteration_num);     % 众数滤波，和medfilt2一样
    % 填充空隙
    bwfill = imfill(bw_m, 'holes');
    
%     清除面积过小的前景（这一段已经在后面由清除边际长度过短的前景代替，使用面积并非一个好的选择）
%     STATS=regionprops(bwfill,'all');
%      for ind_i=1:length(STATS)
%          if(STATS(ind_i).Area<200)     %清除面积小于200的点  需要调整
%              tmp=STATS(ind_i).PixelList;
%              for j=1:numel(tmp)/2
%                  bwfill(tmp(j,2),tmp(j,1))=0;
%              end
%          end
%      end

    %% 5.提取边缘
    edgeim = bwperim(bwfill); %figure;imshow(edgeim);   
    edgeim = filledgegaps(edgeim, 3);
    
    %% 6. 需要剔除毛刺才能交给edgelink，必须保证轮廓为单层像素 (这段代码在细胞贴边时会出错）
    % 加入padarray后已修复 2015.4.28
    % 调用内部函数
    edgeim = Remove_Spurs( edgeim );
    % 除去较小的轮廓
    [LL, nn] = bwlabel(edgeim,8);
    for ii=1:nn
        if sum(sum(LL==ii))<=remove_small
            edgeim(LL==ii) = 0;
        end
    end
    % 保存edgeim图片
    lunkuo = edgeim; % 将其输出到外部空间中保存
    
    %% 7. 将轮廓按逆时针方向连接起来（实际上，由于一些奇怪的形状，连接方式有可能出现错乱，出现顺时针或8字型连接）
    % 由于这个edgelink函数的复杂性，难以修复这个bug，因此选择在后面的操作中间接弥补
    % 这个地方的10像素目前还没被用到，因为有些bug，所以需要自己去除过短的边际
    [rj, cj, ~, ~] = findendsjunctions(edgeim);
    if ~isempty(rj)
        if iteration_num==5
            for kk=1:numel(rj)
                edgeim(rj(kk),cj(kk)) = 0;
            end
        else
            lunkuo = 'error';
            e = [];
            return;
        end
    end
    [edgelist] = CX_edgelink( edgeim );
    % 其实matlab自带函数 bwtraceboundary 就可以实现这个功能！！2015.12.23
%     [edgelist, ~] = edgelink( edgeim, 10 );  % edgelist是每个前景的边缘像素集合，labelededgeim是像素点对应的前景编号
    
    %% 8. 去除过小的边际（这段可作为cell删减的模板）（利用双wihle循环删除cell，实现变长度循环）   
    % 调用内部函数 2016.3.9将这段删去，放到6.进行删除操作
%     edgelist = Remove_Small_Edgelist( edgelist, remove_small); % 去除边际小于 remove_small 的前景

    %% 9.对每个前景轮廓进行分段线性拟合
    seglist = lineseg( edgelist, tolerance );    %%使用线性分段    
    figure;
    imshow(edgeim); % 后续的绘图都是在edgeim前景轮廓上进行操作
%     imshow(zeros(size(edgeim)));
    hold on;
    
    %% 第一部分内容总结
%##########################################################################  

 % 上面9个步骤为第一过程，包含图片读入、预处理、提取轮廓、修正轮廓、分段线性拟合
 % 接下来进入第二过程，即椭圆拟合部分
 % 对每个前景，按照其凹点数目的多少，拟合为 1 个或多个椭圆，这也是这个脚本最终要实现的功能
 
%##########################################################################

%% ========================= 开始对每个前景的循环  ========================= %%
% 定义变量
n = length(seglist);% n为前景数目
% theta = cell(n, 1);    % 角度信息（实际上没用用到）
r = {};        % 轮廓段
e = {};        % 椭圆

for i_f= 1 : n % 前景编号按从上到下扫描的方式
    
%         disp('-------------------');
%         disp(['第',num2str(frame),'帧 拟合前景',num2str(i_f),'...']);
        % 使用短字符临时变量代替 seglist 和 edgelist
        % ------------------ %
        ed = edgelist{i_f};
        if isequal(ed(end,:), ed(1,:)) % 通常ed的最后一位和第一位相同
        	ed = ed(1:end-1,:); % 需要把最后一位去掉，防止重复
        end
        % ------------------ %
        seg = seglist{i_f}; % seg中的最后一个元素=第一个
        if ~isequal(seg(end,:), seg(1,:))
        	seg = [seg; seg(1,:)]; % 如果不等，则需补上一位
        end
       %% 1.分析凹点个数是否正确（使用自己写的edgelink后，可以在之前就判断是否有错误！2015.11.5）
        max_aodian = 3; % 允许连续出现的最多凹点数（要根据数据集的复杂程度决定）
        [ seg, ed, flag_aodian, hasError ] = Analyse_aodian( seg, ed );
%         if hasError
%             lunkuo = 'error';
%             close('1');
%             return;
%         end
 
       %% 4.以下进行椭圆拟合部分
        na = sum(flag_aodian);   % na为凹点总数目 na = ‘num_aodian’
        
       %% 4.1 点数为1或0则拟合全部
        if na<2        
            r{i_f,1} = ed;
        end
        
       %% 4.2 凹点数大于1的情况下则拟合多段
        if na>=2
            % 注意每个变量都要初始化，否则循环时候容易出错
            % xy_axis 是凹点的坐标
            xy_axis=seg(flag_aodian,:);  
            % i_in_e是凹点在edgelist中的序号
            i_in_e = zeros(na,1);      
            for jj=1:na
                i_in_e(jj) = find(ismember(ed ,xy_axis(jj,:),'rows'));
            end
            % 去除过近的凹点，即压制周边7个像素(不去掉也没关系 2015.10.8)
%             [ ind_in_edgelist ] = CX_sup_aodian(numel(ed)/2, ind_in_edgelist, 7);
            % 更新凹点数目
            na = numel(i_in_e); 

%             % 绘制更新后的凹点,加一个方框标记 （可删去）
%             x00 = ed(i_in_e,1);
%             y00 = ed(i_in_e,2);
%             plot(y00, x00, 'bs', 'LineWidth', 1.5);

            % 多边形按凹点分段，然而最后采用原图点进行拟合
            for j2=1:na-1
                r{i_f,j2} = ed(i_in_e(j2):i_in_e(j2+1),:);    
            end
            % 前na段为单段拟合
            r{i_f,na} = [ed(i_in_e(na):end,:);ed(1:i_in_e(1),:)];  

            % 多段拟合调用内部函数
            % 输入该行cell即可，不需要将整个r输入
            [ r flag_zuhe ] = Assign_PixelList( r, i_f, na ); % 将r填入pixellist
            
        end % end if na>=2
        
       %% 4.3 使用轮廓来拟合椭圆
        [~, width]=size(r(i_f,:));
        % 计算每行的cell数目，即前景的拟合数目
        while isempty(r{i_f,width}) 
                width=width-1;
        end     
            
        for k=1:width
            % 拟合方法（较好）调用函数 Prasad_ElliFit
            [~, a, b, x0, y0, alpha, status] = Prasad_ElliFit(r{i_f,k}(:,2), r{i_f,k}(:,1));
            if status && a<=b*5 && size(r{i_f,k},1)>=8;    % 画图的约束条件，不绘制太奇怪的椭圆
                c=cos(alpha);
                s=sin(alpha);
                polar_angle=linspace(0,360,361);
                xq= a*cosd(polar_angle);
                yq= b*sind(polar_angle);
                xn=xq*c-yq*s+x0;
                yn=xq*s+yq*c+y0;
                plot(xn,yn,'r','LineWidth',1.5);
                % 保存椭圆数据,无论椭圆什么样，都保存（后期可修改为只保存status不为0的）
                e{i_f,k}.a=a;
                e{i_f,k}.b=b;
                e{i_f,k}.x0=x0;
                e{i_f,k}.y0=y0;
                e{i_f,k}.alpha=alpha/pi*180; % 图片的y轴是向下的，因此角度逆时针为-，在正常坐标系里逆时针为+
                e{i_f,k}.status=status;
                e{i_f,k}.similar_delete = 0; % 定义相似性，为1则需删除（后续）
                e{i_f,k}.hd = 0;
                e{i_f,k}.num_hypoth = width;
                e{i_f,k}.ind_region = i_f;
                % 保存组合情况表
                if width>1
                    e{i_f,k}.flag_combine = flag_zuhe(k,:);
                else
                    e{i_f,k}.flag_combine = 1;
                end
            else % 不满足条件则status为0
                e{i_f,k}.status = 0;
            end

        end % 一个前景计算结束
        
       %% 4.4 计算豪斯多夫距离
%        disp(['产生假说',num2str(width),'个']);
%        disp('Calculate Hausdorff distance...')
       e(i_f,:) = CX_Hausdorff_Distance(e(i_f,:), r(i_f,:));  % 计算豪斯多夫距离等
%        disp('完成！');

end
hold off;

toc;

end



% ==================== 以下部分为 CX_FIT 内部调用的函数 =================== %
%
% 有些函数原先在外部，后来移入，方便管理
% 全部加入path之后还是放在外面比较整齐 2015.10.7
% ======================================================================= %





















