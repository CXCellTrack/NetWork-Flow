% ����2��ȥ��ë��
function edgeim = Remove_Spurs( edgeim )

flag=1;
edgeim = padarray(edgeim, [1 1], 0);  % ��չ��Ե
maoci_33 = CX_maocijuzhen(0); % ͨ��ë�̾������ö�ٿ��ܲ�����ë��

% tic;
while flag
    pre_edgeim = edgeim;
    STATS=regionprops(pre_edgeim, 'PixelList');    % regionprops������pixellist��������ʱ�����У������Ҫ��edgelink����
    flag=0;
    for ind_i=1:numel(STATS)
        pe=STATS(ind_i).PixelList;      % �����굼��edgelist
        for ind_j=1:size(pe,1)
            px=pe(ind_j,1);
            py=pe(ind_j,2);
            num_33=pre_edgeim(py-1:py+1,px-1:px+1);  % 3*3��Ѱ�Ұ׵㣬��ĿΪ2˵���ǶԽ�ë��

            % -------------------------------
            % �·���������cellfun�жϣ�1.04�룬������(cellfun�����ܲ�����for�ã�)
%                 compare_33 = cellfun(@(x) isequal(num_33, x), maoci_33, 'un', 1);
%                 flag = sum(compare_33);
%                 pre_edgeim(py,px) = ~flag;
            % -------------------------------
            % ԭʼ����������forѭ���жϣ�0.17��
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