function [ seg, ed, flag_aodian, hasError ] = Analyse_aodian( seg, ed, max_aodian )

hasError = 0;

n_seg = size(seg, 1);
% ���Ʒֶ�ͼ��ֻ��Ϊ��debug���㣬����ʵ���ô���
for j=1:n_seg-1                                    % x��������y��������          
    x_j = seg(j,2); y_j = seg(j,1);
    x_k = seg(j+1,2); y_k=seg(j+1,1);
    % ���ߣ��ⲿ�ֿ��Բ��ã�����replot�Ͳ�������
    plot([x_j x_k],[y_j y_k], 'g', 'LineWidth', 1.5); 
end

%% 1.�ж����Էֶε�ÿ���˵㰼�㻹��͹�㣨�ȼ���segΪ��ʱ�����У�

%     1 <--- end-1
%    / 
%   < 
%  2 ---> 3 ...
%
% ��ͼΪʾ��ͼ����123Ϊ͹��ʱ det([v1;v2])>0

flag_aodian = false(1, n_seg);          % Ԥ�����жϰ�����ڵı�ʶ��
for j=1:n_seg-1    
    if j==1 % ��һ������Ҫ���⴦��
        v1=seg(j,:)-seg(n_seg-1,:);       % ǰһ�㵽��ǰ������
        v2=seg(j+1,:)-seg(j,:);      % ��ǰ�㵽��һ������       
    else % ������
        v1=seg(j,:)-seg(j-1,:);     
        v2=seg(j+1,:)-seg(j,:);     % edgelist ����ʱ������ô����
    end
    % ��˺�����������Ĵ�С�ͷ�����Ϣ
    v_multi = det([v1;v2]);                 

    if v_multi<0    % ˵��Ϊ����
        flag_aodian(j)= 1;    % �����¼
    end              
end

%% 2.�жϰ�����Ŀ��������˵���д���ֻ�������ְ취���޸� edgelink ������bug��

reach_max_aodian = 1; % ȷ���Ƿ������������ı�ʾ��
for j=1:n_seg-(max_aodian-1)
    % ��������4�����㣬˵�� edgelist ���ǰ�������ʱ��������
    for ia=0:max_aodian-1
        if flag_aodian(j+ia)==0 % �����������У�ֻҪ��һ���ǰ�������Ϊ��������������Ĵ������
            reach_max_aodian = 0;
        end
    end
    
    if reach_max_aodian
        seg = fliplr(seg')';  % �����з�ת����Ϊ˳ʱ��
        ed = fliplr(ed')';
        flag_aodian = [ ~flag_aodian(1), ~fliplr(flag_aodian(2:end)) ];
        break;
    end
end

for j=1:n_seg-4
    % ��ת������������4�����㣬˵��������8����
    if flag_aodian(j) && flag_aodian(j+1) && flag_aodian(j+2) && flag_aodian(j+3) && flag_aodian(j+4)
        hasError = 1;
        return;
    end
end
        
% �������2����鶼������˵���޴���
% ------------------------- ��ɾ ------------------------------- %
for j=1:numel(flag_aodian)
    if flag_aodian(j)
        plot(seg(j,2), seg(j,1), 'bs', 'LineWidth', 1.5); % ���㻭һ��X����ǣ����Բ�����
    end
end
% --------------------------------------------------------------- %


% %% 1.��� edgelist �Ƿ�Ϊ��ʱ�����У�������ǣ�����ת��Ϊ��ʱ��
% 
% %     1 ---> end-1
% %    / 
% %   < 
% %  2 ---> 3 ...
% %
% % ��ͼΪʾ��ͼ����ʱ������ʱ det([p21;pn1])>0
% p21 = seg(2,:)-seg(1,:);
% pn1 = seg(end-1,:)-seg(1,:);
% % det����������
% if det([p21;pn1])<0     % С��0˵��seg��˳ʱ��ģ���Ҫ��һ��,ͬʱedgelistҲ��Ҫ��ת
%     seg=fliplr(seg')';  % (2015.4.29)������ǰ������ʱedgelink����������bug
%     ed=fliplr(ed')';
% end
% % �������������� edgelist �� seglist ���������ʱ��











