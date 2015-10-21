function  CX_RePlot_Ellipse( dataset )

% ----------------------------------------------------------------------- %
% ��������������»����Ż�����Բ��ȥ���˷ֶ�ֱ�ߺ��Ż�ʱȥ������Բ
% �����ڡ������ͼ������ļ�������������Բͼ
% ���ƿ��ӻ����ٵ�ͼ�����������Բͼ�Ļ�����
% ----------------------------------------------------------------------- %

% if 1
%     dataset = 'competition'; % ѡ��ѵ�����ǲ���
% else
%     dataset = 'training';
% end
[ segpath trackpath ] = getpath( dataset );
% ������Բ����

load([trackpath, '\Pair\Pre_data_new.mat'], 'Ellipse');
lunkuo_addr = [segpath, '\FOI��ȡ����\'];  % ֻ��Ҫ�޸Ĵ˴�
lunkuo_dir = dir([ lunkuo_addr, '*.tif' ]);
output_addr = [trackpath, '\�����ͼ\'];

frame = numel(Ellipse);

for t=1:frame
    disp(['�����',num2str(t),'֡...']);
    pic_name = [ lunkuo_addr, lunkuo_dir(t).name ];
    edgeim = imread(pic_name);
    % ��������ͼƬ��Ϊ�װ�
    figure;
    imshow(edgeim);hold;
    % ������Բ
    for j=1:numel(Ellipse{t})
        e = Ellipse{t}{j};
        alpha1 = e.alpha;
        a = e.a;
        b = e.b;
        x0 = e.x0;
        y0 = e.y0;

        c = cosd(alpha1);
        s = sind(alpha1);
        polar_angle = linspace(0,360,361);
        xq = a*cosd(polar_angle);
        yq = b*sind(polar_angle);
        xn = xq*c-yq*s+x0;
        yn = xq*s+yq*c+y0;
        plot(xn,yn,'g','LineWidth',1.5,'displayname',num2str(j));
        
    end
    hold;
    savename= strcat( output_addr, lunkuo_dir(t).name(1:end-4), '_fit.fig');
    saveas(1,savename);
    close('1');
end

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    