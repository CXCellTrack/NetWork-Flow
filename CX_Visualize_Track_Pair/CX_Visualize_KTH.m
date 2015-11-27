

clear;close all

kthpath = 'E:\datasets\second_edition\competition_datasets\Fluo-N2DH-SIM+\01_KTH_RES\';
kthdir = dir([kthpath,'*.tif']);
savepath = 'E:\datasets\second_edition\competition_datasets\Fluo-N2DH-SIM+\01_KTH_RES_TRACK\';

%% ����track.txt
KTH_track = load([ kthpath, 'res_track.txt']); % ����KTH��track.txt
KTH_track(:,2:3) = KTH_track(:,2:3) + 1;

% ͳ���¼�
divide = KTH_track(KTH_track(:,4)~=0,:);
other = KTH_track(KTH_track(:,4)==0,:);

%% ѡ����ɫ
color = colormap(hsv);close('1'); 
color_numbel = numel(color)/3;
rng(1);
color = color(randperm(color_numbel),:);
nc = size(color,1);
inten2color = zeros(500,1); % �Ҷ�ӳ�䵽��ɫ������


for n=1:numel(kthdir)
    
    disp(['  ����ͼƬ ',num2str(n),'...']);
    im = imread([kthpath, kthdir(n).name]);
    imd1 = zeros(size(im));imd2 = imd1;imd3 = imd1;
    % ͳ��im�е�����  
    stats = regionprops(im, 'Centroid');
    % ��һ֡������ɫ
    if n==1 && numel(stats)<=nc
    	inten2color(1:numel(stats)) = 1:numel(stats); 
    end
    
    %% ������ɫ
    for j=1:numel(stats)
        if inten2color(j)==0 % ˵�����³��ֵ�
            if j<=nc % û������������Ӧ
                inten2color(j) = j;
            else % �����������������һ��
                inten2color(j) = randi(nc);
            end
        end
        imd1(im==j) = color(inten2color(j),1);
        imd2(im==j) = color(inten2color(j),2);
        imd3(im==j) = color(inten2color(j),3); 
        imd = cat(3, imd1,imd2,imd3);
    end
    imshow(imd); hold on;

    %% ͳ�Ʒ��Ѳ����
    ind = find(divide(:,2)-1==n);
    if ~isempty(ind)
        assert(mod(numel(ind),2)==0); % ����Ϊż����
        for j=ind'
            son = divide(j,1);
            father = divide(j,4);
%             inten2color(son) = inten2color(father); % �̳���ɫ(���׿������)
            inten2color(son) = randi(nc); % ���⿪һ����ɫ
            % ���б��
            e.x0 = stats(father).Centroid(1);
            e.y0 = stats(father).Centroid(2);
            CX_plot_event( e, 'd');
        end
    end
    
    %% ͳ�Ƴ��ֲ����
    if n~=1
        ind = find(KTH_track(:,2)==n);
        if ~isempty(ind)
            for j=ind'
                tmpe = KTH_track(j,1);
                e.x0 = stats(tmpe).Centroid(1);
                e.y0 = stats(tmpe).Centroid(2);
                
                if any(divide(:,4)==KTH_track(j,4)) % �������ϸ������
                    CX_plot_event( e, 'dson');
                else
                    CX_plot_event( e, 's'); % ���ֱ��
                end
            end
        end
    end
    
    %% ͳ����ʧ�����
    if n~=numel(kthdir)
        ind = find(other(:,3)==n);
        if ~isempty(ind)
            for j=ind'
                tmpe = other(j,1);
                if any(divide(:,4)==tmpe) % �����ĸϸ��������ʧ
                    continue;
                end
                
                e.x0 = stats(tmpe).Centroid(1);
                e.y0 = stats(tmpe).Centroid(2);
                CX_plot_event( e, 't');
            end
        end
    end
    hold off
    
    %% ����
    saveas(1,[savepath, strrep(kthdir(n).name,'tif','fig')]);
    close(1);
    
end
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        

            