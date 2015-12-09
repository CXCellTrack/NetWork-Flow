function superpixels = calculate_sp_info( row, label_zuhe, flag_zuhe, stats, superpixels )

% ��ÿ����˵������Ϣ���൱����Բ�е���Ϣ��

for j=1:numel(label_zuhe)

    % ----------------- ��������pixellist ----------------- % 
    sum_pixellist = [];
    for k=label_zuhe{j}
        sum_pixellist = [ sum_pixellist; stats(k).PixelList];
    end
    % ------------------------------------------------------ %
    % �ȷ�5��ѡ��
    superpixels{row}{j,1}.label = label_zuhe{j};
    superpixels{row}{j,1}.pixellist = sum_pixellist;
    superpixels{row}{j,1}.centroid = [mean(sum_pixellist(:,1)), mean(sum_pixellist(:,2))];
    superpixels{row}{j,1}.flag_combine = flag_zuhe(j,:);
    superpixels{row}{j,1}.ind_region = row;
    superpixels{row}{j,1}.num_hypoth = numel(label_zuhe);
end