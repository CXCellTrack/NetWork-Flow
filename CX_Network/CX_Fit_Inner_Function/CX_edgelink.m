function edgelist = CX_edgelink( edgeim )

% ���Լ�д��edgelink�����������ҵ��Ǹ�
% ����ٶ���΢��һ��
% ���������� findjunction ���ж��Ƿ���ʮ�ֽ���ĸ������
% 2015.11.5


% profile on
% tic
% im = imread('E:\datasets\second_edition\competition_datasets\Fluo-N2DH-SIM+\01_4-16_seg\FOI��ȡ����\t001.tif');
% [edgelist, ~] = edgelink( im, 10 );toc

[ROWS, COLS] = size(edgeim);
stats = regionprops(edgeim, 'PixelList');
edgelist = cell(size(stats));

roff = [-1 -1 -1  0  0  0  1  1  1]; % 3*3�����ڵ�����ƫ��
coff = [-1  0  1 -1  0  1 -1  0  1];
alloff = abs(roff) + abs(coff); % ƫ��֮�ͣ������ж���ʮ�����ӻ��ǶԽ�����

for n=1:numel(stats)
    xyco = stats(n).PixelList*[0 1;1 0];
    center = xyco(1,:); % ��������
    edgelist{n} = [ edgelist{n}; center]; % ����켣
    
    while 1
        % 1 2 3
        % 4 5 6
        % 7 8 9
        r = center(1) + roff;
        c = center(2) + coff;
        ind = find((r>=1 & r<=ROWS) & (c>=1 & c<=COLS));

        
        if size(edgelist{n},1)==1 % ��һ�������Һ�����ȥ�ң�ֻ������8��9ȥ��
            if any(ind==8) && edgeim(r(8),c(8))==1
                edgelist{n} = [ edgelist{n}; r(8),c(8) ];
                center = [r(8),c(8)];
                direction = 10 - 8; % ��ǰ���Ķ����¸�������˵�ķ���
            elseif any(ind==9) && edgeim(r(9),c(9))==1
                edgelist{n} = [ edgelist{n}; r(9),c(9) ];
                center = [r(9),c(9)];
                direction = 10 - 9;
            end
        else
            % ���ǵ�һ���㣬����������Ѱ���¸���
            
            tmp = [];
            switch direction % fromΪ��·���ϵģ�����ѡ
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
            % ��ind��ɾȥfrom�е�  
            % ʹ��setdiff0.4�룬ʹ��mysetdiff0.2�룬������
%             ind = mysetdiff(ind, from);
            for ii=ind
                if edgeim(r(ii),c(ii))==1 && ii~=5 && all(from~=ii)
                    tmp = [tmp, ii];
                end
            end
            if numel(tmp)==2 % �ҵ�2���µ㣬˵��һ��ʮ�֣�һ���Խǣ��Ǿ�ѡ��ʮ��
                ti = find(alloff(tmp)==1);
            else
                ti = 1; % ֻ��һ�����û��ѡ
            end
            ii = tmp(ti);
            edgelist{n} = [ edgelist{n}; r(ii),c(ii) ];
            center = [r(ii),c(ii)];
            direction = 10 - ii;
        end
        
        % �����һ������Ϊ��һ���㣬˵���ƻ�����
        if isequal(center, edgelist{n}(1,:))
            break;
        end
        
    end
    
end

                    
% toc
% profile viewer     
            
            
            
                    
                    
                    
                
    
    
    
    
    
    
    
    