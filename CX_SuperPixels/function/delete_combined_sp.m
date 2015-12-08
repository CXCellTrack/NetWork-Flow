function [ label_zuhe flag_zuhe ] = delete_combined_sp( fore, label_nearby, label_zuhe, flag_zuhe )

% sp ���� superpixels
% fore��ǰ���е�����basic sp
% label_zuhe��ȫ��ϣ���2^n-1��
% label_nearby��ÿ��label���ڵ�label����

% ����ɾ�������ڵ���Ϸ�ʽ
% ======================================= %

% �������Щ�ǲ����ڵģ�Ҫȥ��
for bsp=fore % ��fore�е�basic_superpixel����
    for i_csp=numel(fore)+1:numel(label_zuhe) % �����е�combine_superpixel����
        csp = label_zuhe{i_csp};
        if any(csp==bsp) % �������а���bsp
            rest = mysetdiff(csp,bsp); % ����bspʣ�µĵ���bsp������������
            % A������ʹ�ü��������ٶ���΢��һ��
            % ����ཻΪ�գ�˵��rest��bsp�������ڣ���ɾ�����ּ�˵
%                 if isempty( intersect(rest, label_nearby{bsp}) ) 
%                     label_zuhe{i_csp} = [];
%                 end
            % B������arraufun����Ƚϣ��ٶȿ�
            tmp = arrayfun(@(x) any(rest==x),label_nearby{bsp},'un',1);
            if all(tmp==0)
                label_zuhe{i_csp} = [];
            end
        end
    end
end

ind_save = ~isemptycell(label_zuhe);
label_zuhe = label_zuhe(ind_save);
flag_zuhe = flag_zuhe(ind_save,:);

