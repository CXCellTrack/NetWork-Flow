%######################################
%
% 2015.5.30 CX on desk
% ���ã�����������ڽ�������ɵ�ϸ����ɫ�����ڹ۲�
% ���ݴ洢������ɫ��ͼƬ�����ڡ����ӻ����ٱ�ǡ�
% ������ϵ������ CX_RePlot_Ellipse ���Ƹ��º����Բ
%          ���� CX_fill_color �� CX_fill_color_event ���и��ٱ��
%###################################### 

%% ���ƿ��ӻ�����
clear;close all;

% ָ�����ĸ����ݼ��Ͻ��м��㣨train or test��
if 1
    dataset = 'competition';
else
    dataset = 'training';
end
[ ~, trackpath ] = getpath( dataset );

% ����ȫ�ֱ���
global Fij Fit Fid Fiv Fmj Fsj;
global conflict_fij conflict_pair_last_xy conflict_pair_next_xy n;

load([trackpath, '\Pair\Pre_data_New.mat']);
fig_addr = [trackpath, '\�����ͼ\'];
 
if 0
    disp('  ������ʵ���̱�������...');
    load([trackpath, '\GT\GT_after_hand_tune\GT_Flow_Variables_New.mat']);
    output_addr = [trackpath, '\GT\GT_after_hand_tune\'];
    s_frame = 1;
    e_frame = numel(Fmj);
else
    disp('  ���뾭��SSVM learing��õ��ķ��䷽��');
    track_data_addr = [trackpath, '\�ṹ��ѧϰ\asl_track.mat'];
    output_addr = [trackpath, '\�ṹ��ѧϰ\asl_track\'];
    mkdir(output_addr)
    load( track_data_addr );
    % ������������̱������� s_frame �� e_frame
    s_frame = sum(cellfun(@isempty, Fsj));
    e_frame = numel(Fsj);
end
frame = e_frame - s_frame + 1;

color = colormap(hsv);close('1'); % ѡ����ɫ
nc = size(color,1);
rng(1)
color = color(randperm(nc),:);

%% �����ͼ�װ壬�������ͼ
fig_dir = dir([ fig_addr, '*.fig' ]);
% ���ñ���ͼƬ���ļ��У���������е�figͼƬ
delete([output_addr,'*.fig']);
% ���������ͼ�ļ����ڵ�ͼƬ������Ŀǰ��frameСʱ����Ҫ�ػ����ͼƬ
if 0
    CX_RePlot_Ellipse( dataset );
end
% ��Բ�Ƿ���Ҫ�߿��ޱ߿��Ư�����ʺ���ͼ���б߿����ڲ��debug
need_contour = true;

%% ��ʼ����ѭ������
for t = s_frame:e_frame
    fig_name = [ fig_addr, fig_dir(t).name ];
    if need_contour
        openfig(fig_name, 'new'); %'invisible'
    else
        sizeim = [ ]; % �ޱ߿�ģʽ����Ҫ��дͼƬ��С
        im = zeros(sizeim);
        imshow(im);
    end
    
    disp(['  ���ڴ���', fig_name, '...']);
    hold on 
    %##### ��һ֡��Ҫ���������ɫ #####
    if t==s_frame
        for j=1:n(t)
            % ������
            [ eventIn eventOut ] = CX_CheckInOut( t, j );
            if sum(eventOut)>0
                Ellipse{t}{j}.color_index = randi(nc); % �г��ڵ�ϸ���ͷ�����ɫ
            end
        end
    end
        
    for j=1:n(t)
        %#################### ��ɫ #######################
        if isfield( Ellipse{t}{j}, 'color_index' )
            % �����Բ��2����ɫ����˵����Ϊmerge��Ǩ�����ģ�������ɫ
            % ע�⣬���ﲻ�ܽ�merge������ɫ����Ϊ���е�����ʱ��û�ж��Ƿ�Ϊmerge
            % merge��ɫ�Ĳ������ж�merge���������Ѿ�������������ɫ�ˣ�
            assert(numel( Ellipse{t}{j}.color_index )<=2);
            if numel( Ellipse{t}{j}.color_index )==2
                CX_fill_color( Ellipse{t}{j}, [1 1 1] ); % ���merge֮����move�����
            elseif numel( Ellipse{t}{j}.color_index )==1
                CX_fill_color( Ellipse{t}{j}, color, Ellipse{t}{j}.color_index );
            end
        end
    
        %################# ���ݲ�ͬ���¼�����ɫ���ݸ���� ##################
        % �����ڳ���
        [ eventIn eventOut ] = CX_CheckInOut( t, j );
        
        %########## ����¼��������ڵ�һ֡ ##########
        if t~=s_frame
            ev = find(eventIn);
            if ~isempty(ev) % �ӿ��ٶ�
                switch ev
                    case 3 % ��ϸ��
                        CX_plot_event( Ellipse{t}{j}, 'dson' );
                    case 4 % ���ѵõ���Сϸ��
                        CX_plot_event( Ellipse{t}{j}, 'vson' );
                        
                    case 5 % merge
                        mm = find(Fmj{t}(j,:)==1);
                        source = candidate_k_last{t}{j,mm};
                        Ellipse{t}{j}.color_index = [ Ellipse{t-1}{source(1)}.color_index, Ellipse{t-1}{source(2)}.color_index ];
                        CX_plot_event( Ellipse{t}{j}, 'm');   % �ںϵ�����ϸ�����
                        CX_fill_color( Ellipse{t}{j}, [1 1 1] );  % �ںϵõ�����Ҫ����ɫ

                    case 6 % appear
                        CX_plot_event( Ellipse{t}{j}, 's');  % �³���ϸ�����
                        Ellipse{t}{j}.color_index = randi(nc); 
                        CX_fill_color( Ellipse{t}{j}, color, Ellipse{t}{j}.color_index ); % ���ֱȽ����⣬�³��ֵ���Ҫ��������ɫ
                end
            end
        end
        
        %########## �����¼������������һ֡ ##########
        if t~=e_frame
            ev = find(eventOut);
            if ~isempty(ev) % �ӿ��ٶ�
                switch ev
                    case 1 % move
                        mm = find(Fij{t}(j,:)==1);
                        e_next = candidate_fij{t}(j,mm);
                        Ellipse{t+1}{e_next}.color_index = Ellipse{t}{j}.color_index;  

                    case 2 % disappear
                        CX_plot_event( Ellipse{t}{j}, 't' );

                    case 3 % divide
                        mm=find(Fid{t}(j,:)==1);
                        % �ҵ�2����ϸ���������Ǽ̳���ɫ�����ñ�ǣ���ɫ�����µģ�
                        son = candidate_k_next{t}{j,mm};

                        Ellipse{t+1}{son(1)}.color_index = randi(nc);
                        Ellipse{t+1}{son(2)}.color_index = randi(nc);

                        CX_plot_event( Ellipse{t}{j}, 'd' );

                    case 4 % split
                        mm=find(Fiv{t}(j,:)==1);
                        % �ҵ�2������ϸ���������Ǽ̳���ɫ
                        split = candidate_k_next{t}{j,mm};
                        if numel(Ellipse{t}{j}.color_index)==2 % ����ҵ�2����ɫ�͸�����ȡ������˵����merge��split�ˣ��д���
                            Ellipse{t+1}{split(1)}.color_index = Ellipse{t}{j}.color_index(1);
                            Ellipse{t+1}{split(2)}.color_index = Ellipse{t}{j}.color_index(2);
                        else
                             Ellipse{t+1}{split(1)}.color_index = randi(nc);
                             Ellipse{t+1}{split(2)}.color_index = randi(nc); % �������
                        end
                        CX_plot_event( Ellipse{t}{j}, 'v' );  
                        
                    case 5 % merge
                        CX_plot_event( Ellipse{t}{j}, 'mson' );
                end
            end
        end
    end
    hold off
    savename = [output_addr, strrep(fig_dir(t).name,'tif','fig')];
    saveas(1,savename);
    close('1');
end

%% ͼƬ��ʾ
if 0
    screen_size = get(0,'ScreenSize');
    colored_fig_dir = dir([ output_addr, '*.fig' ]);  
    for t = 1:20
        fig_name = [ output_addr, colored_fig_dir(t).name ];
        openfig(fig_name, 'new', 'visible');
        set(gcf, 'Position', screen_size);
    end
end




