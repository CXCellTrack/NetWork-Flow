function [ label_zuhe flag_zuhe ] = delete_uncombined_sp( fore, label_nearby, label_zuhe, flag_zuhe, bsp_stats, new_labels )

% sp ���� superpixels
% fore��ǰ���е�����basic sp
% label_zuhe��ȫ��ϣ���2^n-1��
% label_nearby��ÿ��label���ڵ�label����

% ����ɾ�������ڵ���Ϸ�ʽ
% ======================================= %


%% A������ʹ��ͼ���ж�����Ƿ�������ٶ�������0.059�룩������֤��ȷ��
% tic
for i_csp=1:numel(label_zuhe) % ��fore�е�����superpixel����
    % ��ǰsp
    csp = label_zuhe{i_csp}; 
    % ��ͼ�װ�
    bb = zeros(size(new_labels)); 
    for bsp=csp
        pl = bsp_stats(bsp).PixelList;
        for i_xy=1:size(pl,1) % ��bsp����ͼ�ϣ��ټ��ͼ����ͨ��
            bb(pl(i_xy,1), pl(i_xy,2)) = 1;
        end
    end
    [~, nlabel] = bwlabel(bb, 4);
    
    if numel(csp)==1 % ����bsp��������ֻ��һ�飬�������ǰ�������
        assert(nlabel==1)
    end
    if nlabel~=1
        label_zuhe{i_csp} = [];
    end
end

ind_save = ~isemptycell(label_zuhe);
label_zuhe = label_zuhe(ind_save);
flag_zuhe = flag_zuhe(ind_save,:);
% toc

if 0
%% B����������Ҫ��ͼ��ʹ�����ڹ�ϵ���ж�����Ƿ�������ٶȿ죬����ʱ�����⣩% 2015.12.10 ��������һ������0.033�룩
tic
% �������Щ�ǲ����ڵģ�Ҫȥ��
for bsp=fore % ��fore�е�basic_superpixel����
    % bsp���������򹹳ɵļ���
    S_nb = [label_nearby{bsp}, bsp]; 
    for i_csp=numel(fore)+1:numel(label_zuhe) % �����е�combine_superpixel����(basic_sp����Ҫ�����ͨ��)
        csp = label_zuhe{i_csp};
        if any(csp==bsp) % �������а���bsp
            neighbor_in_csp = intersect(S_nb, csp); % csp�а�����bsp�������ھ�
            rest = mysetdiff(csp, S_nb); % ʣ�µĲ���
            if isempty(rest) % ���ûʣ�£�˵�����У�ֱ������
                continue;
            end
            % ���S_nb�����е�bsp����rest����������˵��restΪ����һ�飬���ʧ��
            delete_flag = true;
            for kk=neighbor_in_csp
                tmp = arrayfun(@(x) any(rest==x),label_nearby{kk},'un',1);
                if any(tmp==1) % ֻҪbsp�����ھ���һ�������ӵ�rest����˵����ͨ
                    delete_flag = false;
                    break;
                end
            end
            if delete_flag
                label_zuhe{i_csp} = [];
            end
        end
    end
end

ind_save = ~isemptycell(label_zuhe);
label_zuhe = label_zuhe(ind_save);
flag_zuhe = flag_zuhe(ind_save,:);
toc

end
