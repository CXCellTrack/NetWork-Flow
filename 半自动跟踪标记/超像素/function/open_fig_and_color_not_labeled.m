function open_fig_and_color_not_labeled( label2e, SuperPixel, fig_dir, t )

[ ~, trackpath ] = getpath( 'training' );
% ������ͼ�Ĵ�С
screen_size = get(0,'ScreenSize');
first_pos = [1, 50, screen_size(3)/2, screen_size(4)-130];
second_pos = [screen_size(3)/2, 50, screen_size(3)/2, screen_size(4)-130];

% ����ͼƬ
openfig( [ trackpath, '\GT\label_and_e\', fig_dir(t).name] );
if gcf==1
    set(gcf, 'Position', first_pos);
elseif gcf==2
    set(gcf, 'Position', second_pos);
end

set(gca, 'position',[0 0 1 1]);
h = get(gcf, 'children');
% ------------------------------------------------------------------- %
% ��û���Զ�����ϵ���Բ����Ϊ��ɫ��*�����Ϊ��ɫ���Լ�ǿ�Ӿ�����
% 1��������ɫ��
label_not_attached = find(label2e{t,1}==0);
for lab=label_not_attached'
    h_label = findobj(h, 'marker','*','DisplayName',num2str(lab) );
    set(h_label, 'color', 'b');
end
% 2�����ƺ�ɫԲȦ
bsp = cell2mat( cellfun(@(x) x.label, SuperPixel{t},'un',0)' );
maxbsp = max(bsp); % �ҵ����ͼ������basic sp��Ҳ����basic sp������Ŀ��
sp_not_labed = setdiff( 1:maxbsp, label2e{t}(:,2) );
for jj=sp_not_labed
    h_e = findobj(h, 'color','g','DisplayName',num2str(jj) );
    set(h_e, 'color', 'r'); % ��ɫ
end
% ------------------------------------------------------------------- %
