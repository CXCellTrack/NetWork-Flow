function superpixels = calculate_sp_info( row, label_zuhe, flag_zuhe, stats, superpixels )

% 对每个假说计算信息（相当于椭圆中的信息）

for j=1:numel(label_zuhe)

    % ----------------- 求解区域的pixellist ----------------- % 
    sum_pixellist = [];
    for k=label_zuhe{j}
        sum_pixellist = [ sum_pixellist; stats(k).PixelList];
    end
    % ------------------------------------------------------ %
    % 先放5个选项
    superpixels{row}{j,1}.label = label_zuhe{j};
    superpixels{row}{j,1}.pixellist = sum_pixellist;
    superpixels{row}{j,1}.centroid = [mean(sum_pixellist(:,1)), mean(sum_pixellist(:,2))];
    superpixels{row}{j,1}.flag_combine = flag_zuhe(j,:);
    superpixels{row}{j,1}.ind_region = row;
    superpixels{row}{j,1}.num_hypoth = numel(label_zuhe);
end