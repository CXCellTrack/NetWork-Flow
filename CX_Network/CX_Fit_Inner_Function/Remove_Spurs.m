% 函数2：去除毛刺
function edgeim = Remove_Spurs( edgeim )

flag=1;
edgeim = padarray(edgeim, [1 1], 0);  % 扩展边缘
maoci_33 = CX_maocijuzhen(0); % 通过毛刺矩阵表来枚举可能产生的毛刺

% tic;
while flag
    pre_edgeim = edgeim;
    STATS=regionprops(pre_edgeim, 'PixelList');    % regionprops产生的pixellist并不按逆时针排列，因此需要用edgelink函数
    flag=0;
    for ind_i=1:numel(STATS)
        pe=STATS(ind_i).PixelList;      % 将坐标导入edgelist
        for ind_j=1:size(pe,1)
            px=pe(ind_j,1);
            py=pe(ind_j,2);
            num_33=pre_edgeim(py-1:py+1,px-1:px+1);  % 3*3内寻找白点，数目为2说明是对角毛刺

            % -------------------------------
            % 新方法，采用cellfun判断，1.04秒，更加慢(cellfun的性能并不比for好！)
%                 compare_33 = cellfun(@(x) isequal(num_33, x), maoci_33, 'un', 1);
%                 flag = sum(compare_33);
%                 pre_edgeim(py,px) = ~flag;
            % -------------------------------
            % 原始方法，采用for循环判断，0.17秒
            for ind_k=1:numel(maoci_33)
                if isequal(num_33, maoci_33{ind_k})
                    flag = 1;
                    pre_edgeim(py,px) = 0;
                end
            end
            % -------------------------------
        end
    end
    edgeim = pre_edgeim;
end
edgeim = edgeim(2:end-1,2:end-1);
% toc;imshow(edgeim);

end