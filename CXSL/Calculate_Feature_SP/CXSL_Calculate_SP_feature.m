function  SuperPixel = CXSL_Calculate_SP_feature( SuperPixel, segpath, n )
% structured learing���֡���������ȡ
% ��ȡ����Superpixel������������ 
%

% ��Ҫ����ԭʼͼƬ�ĵ�ַ
last = max(strfind(segpath, '\'));
raw_addr = segpath(1:last+2);  
raw_dir = dir([ raw_addr, '\*.tif' ]);

frame = numel(SuperPixel);

%% ��ʼѭ��������ÿ��SP������
for t=1:frame % ÿ֡
    if isempty(SuperPixel{t}) % ��ȥ��cell
        continue;
    end
    disp(['  �����',num2str(t),'֡...']);
    im = imread([ raw_addr, '\', num2str(raw_dir(t).name) ]);
    im_db = mat2gray(im);

    for j=1:n(t) % ÿ��SP 
        SP = SuperPixel{t}{j};
        % baseboard ��ͼ�װ壬�ڵװ��ϻ���SP��λ��
        bb = zeros(size(im));
        for pl=1:size(SP.pixellist,1)
            bb(SP.pixellist(pl,2), SP.pixellist(pl,1)) = 1;
        end % imshow(bb)
        % ����SP������
        spfeature = regionprops(bb, im_db, 'MeanIntensity','PixelValues',... % �Ҷ�����
            'Area','Centroid',... % ���������
            'Eccentricity','MajorAxisLength','MinorAxisLength','Orientation',... % ����Բ����
            'Extent','Solidity');
        spfeature.DeviaIntensity = std(spfeature.PixelValues); % ����ӻҶȱ�׼��
        spfeature.SumIntensity = sum(spfeature.PixelValues); % �ҶȺ�
        spfeature.HistIntensity = hist(spfeature.PixelValues, 16)'/spfeature.Area; % �Ҷ�ֱ��ͼ
        % ���߾���
        spfeature.Dist2border = min([ size(im_db)-[SP.centroid(2), SP.centroid(1)], SP.centroid(2), SP.centroid(1) ]);
        % �Ƴ�����ֵ
        spfeature = rmfield(spfeature,'PixelValues');
        % ��ֵ��ɾ���������ԣ�ֻ��������
        SuperPixel{t}{j} = spfeature;
        
    end
end % end frame

    
            



