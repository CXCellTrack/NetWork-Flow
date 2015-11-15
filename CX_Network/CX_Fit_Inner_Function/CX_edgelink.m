function edgelist = CX_edgelink( edgeim )

% 用自己写的edgelink来代替网上找的那个
% 这个速度稍微快一点
% 不过得先用 findjunction 来判断是否有十字交叉的复杂情况
% 2015.11.5


% profile on
% tic
% im = imread('E:\datasets\second_edition\competition_datasets\Fluo-N2DH-SIM+\01_4-16_seg\FOI提取轮廓\t001.tif');
% [edgelist, ~] = edgelink( im, 10 );toc

[ROWS, COLS] = size(edgeim);
stats = regionprops(edgeim, 'PixelList');
edgelist = cell(size(stats));

roff = [-1 -1 -1  0  0  0  1  1  1]; % 3*3矩阵内的坐标偏差
coff = [-1  0  1 -1  0  1 -1  0  1];
alloff = abs(roff) + abs(coff); % 偏差之和，用来判断是十字连接还是对角连接

for n=1:numel(stats)
    xyco = stats(n).PixelList*[0 1;1 0];
    center = xyco(1,:); % 中心坐标
    edgelist{n} = [ edgelist{n}; center]; % 加入轨迹
    
    while 1
        % 1 2 3
        % 4 5 6
        % 7 8 9
        r = center(1) + roff;
        c = center(2) + coff;
        ind = find((r>=1 & r<=ROWS) & (c>=1 & c<=COLS));

        
        if size(edgelist{n},1)==1 % 第一个点往右和右下去找（只可能往8和9去）
            if any(ind==8) && edgeim(r(8),c(8))==1
                edgelist{n} = [ edgelist{n}; r(8),c(8) ];
                center = [r(8),c(8)];
                direction = 10 - 8; % 当前中心对于下个中心来说的方向
            elseif any(ind==9) && edgeim(r(9),c(9))==1
                edgelist{n} = [ edgelist{n}; r(9),c(9) ];
                center = [r(9),c(9)];
                direction = 10 - 9;
            end
        else
            % 不是第一个点，则在领域中寻找下个点
            
            tmp = [];
            switch direction % from为来路边上的，不能选
                case 1
                    from = [1 2 4];
                case 2
                    from = [1 2 3];
                case 3
                    from = [6 2 3];
                case 4
                    from = [1 4 7];
                case 6
                    from = [3 6 9];
                case 7
                    from = [4 7 8];
                case 8
                    from = [9 7 8];
                case 9
                    from = [9 6 8];
            end
            % 从ind中删去from中的  
            % 使用setdiff0.4秒，使用mysetdiff0.2秒，都很慢
%             ind = mysetdiff(ind, from);
            for ii=ind
                if edgeim(r(ii),c(ii))==1 && ii~=5 && all(from~=ii)
                    tmp = [tmp, ii];
                end
            end
            if numel(tmp)==2 % 找到2个新点，说明一个十字，一个对角，那就选择十字
                ti = find(alloff(tmp)==1);
            else
                ti = 1; % 只有一个点就没得选
            end
            ii = tmp(ti);
            edgelist{n} = [ edgelist{n}; r(ii),c(ii) ];
            center = [r(ii),c(ii)];
            direction = 10 - ii;
        end
        
        % 如果下一个中心为第一个点，说明绕回来了
        if isequal(center, edgelist{n}(1,:))
            break;
        end
        
    end
    
end

                    
% toc
% profile viewer     
            
            
            
                    
                    
                    
                
    
    
    
    
    
    
    
    