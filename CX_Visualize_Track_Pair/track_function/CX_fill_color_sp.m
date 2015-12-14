function  im = CX_fill_color_sp( sp, im, color, index )

% 将sp话上颜色

if ~exist('index', 'var')
    thiscolor = color; % 也可以直接输入color不输入index
else  
    assert(numel(index)==1)
    thiscolor = color(index,:);
end

for ii=1:size(sp.pixellist,1)
    im(sp.pixellist(ii,2),sp.pixellist(ii,1),:) = thiscolor;
end

end


