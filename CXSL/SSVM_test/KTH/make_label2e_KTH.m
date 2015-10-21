function stats = make_label2e_KTH( KTH_RES_PATH, e_used, Ellipse )
%
% ��kth�е�ϸ�����ҵķ����е���Բ��˵һһ��Ӧ�ϣ��Է�����㾫��
%

KTH_RES = dir([KTH_RES_PATH, '*.tif']);
stats = cell(numel(KTH_RES),1);

% ���û������Ҫ����һ��stats��ʱ�仨�ѽϳ���
for t=1:numel(KTH_RES)
    im = imread([KTH_RES_PATH, KTH_RES(t).name]);
    erow = e_used(t,e_used(t,:)~=0); % ��ȡ��ʵ�õ�����Բ��˵

    stats{t} = regionprops(im, 'PixelList');
    for u=1:numel(stats{t})
        stats{t}(u).e = []; % һ����˵�㲻����ͬʱλ��2�����򣬵�1�������ڿ����кü�����˵��
    end

    disp(['  �����',num2str(t),'��ͼƬ...']);
    tic
    for j=1:numel(erow) % ��ǵ�ǰ�����Ӧ����Բ��˵
        tmpe = erow(j);
        xy = round([Ellipse{t}{tmpe}.x0, Ellipse{t}{tmpe}.y0]);
        % ���Ҽ�˵���ĵ��Ƿ���������
        for u=1:numel(stats{t})
            if isempty(stats{t}(u).PixelList) % �˻Ҷ�û����������
                continue;
            end
            flag = ismember(stats{t}(u).PixelList, xy, 'rows');
            if any(flag==1)
                stats{t}(u).e = [stats{t}(u).e; tmpe];
                break;
            end
        end
    end
    toc
%     num_e = arrayfun(@(x)numel(x.e), stats{t}); % ÿ��ǰ���а�������Բ����    
end
    
    
    
    
    