%######################################
%
% 2015.5.30 CX on desk
% ���ã�����������ڽ�������ɵ�ϸ����ɫ�����ڹ۲�
% ���ݴ洢������ɫ��ͼƬ�����ڡ����ӻ����ٱ�ǡ�
% ������ϵ������ CX_RePlot_Ellipse ���Ƹ��º����Բ
%          ���� CX_fill_color_sp �� CX_fill_color_event ���и��ٱ��
%###################################### 

%% ���ƿ��ӻ�����
clear;close all;

% ָ�����ĸ����ݼ��Ͻ��м��㣨train or test��
if 1
    dataset = 'competition';
else
    dataset = 'training';
end
[ trainpath, trackpath ] = getpath( dataset );

% ����ȫ�ֱ���
global Fij Fit Fid Fiv Fmj Fsj;
global conflict_fij conflict_pair_last_xy conflict_pair_next_xy n;

load([trackpath, '\Pair\Pre_data_New.mat']);
 
if 0
    disp('  ������ʵ���̱�������...');
    load([trackpath, '\GT\GT_after_hand_tune\GT_Flow_Variables_New.mat']);
    output_addr = [trackpath, '\GT\GT_after_hand_tune\'];
    s_frame = 1;
    e_frame = numel(Fmj);
else
    disp('  ���뾭��SSVM learing��õ��ķ��䷽��');
    track_data_addr = [trackpath, '\�ṹ��ѧϰ\Tracking_Data.mat'];
    output_addr = [trackpath, '\Pair\���ӻ����ٱ��\'];
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

%% �����ͼ�װ�
png_addr = [trainpath, '\�����ػҶȱ�ǩ\'];
png_dir = dir([ png_addr, '*.png' ]);
% ���ñ���ͼƬ���ļ��У���������е�figͼƬ
delete([output_addr,'*.fig']);

%% ��ʼ����ѭ������
for t = s_frame:e_frame
    pic_name = [ png_addr, png_dir(t).name ];
%     openfig(fig_name, 'new'); %'invisible'
    origin = imread(pic_name);
    im = ones([ size(origin), 3 ]);
    eval_str = '';
    
    disp(['  ���ڴ���', pic_name, '...']);
    %##### ��һ֡��Ҫ���������ɫ #####
    if t==s_frame
        for j=1:n(t)
            % ������
            [ eventIn eventOut ] = CX_CheckInOut( t, j );
            if sum(eventOut)>0
                SuperPixel{t}{j}.color_index = randi(nc); % �г��ڵ�ϸ���ͷ�����ɫ
            end
        end
    end
        
    for j=1:n(t)
        %#################### ��ɫ #######################
        if isfield( SuperPixel{t}{j}, 'color_index' )
            % �����Բ��2����ɫ����˵����Ϊmerge��Ǩ�����ģ�������ɫ
            % ע�⣬���ﲻ�ܽ�merge������ɫ����Ϊ���е�����ʱ��û�ж��Ƿ�Ϊmerge
            % merge��ɫ�Ĳ������ж�merge����������ɫ������Щ��merge�����Ǩ�Ƶģ�
%             assert(numel( SuperPixel{t}{j}.color_index )<=2);
            if numel( SuperPixel{t}{j}.color_index )>=2
                im = CX_fill_color_sp( SuperPixel{t}{j}, im, [1 1 1] ); % ���merge֮����move�����
            elseif numel( SuperPixel{t}{j}.color_index )==1
                im = CX_fill_color_sp( SuperPixel{t}{j}, im, color, SuperPixel{t}{j}.color_index );
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
                        eval_str = CX_plot_event_sp( eval_str, SuperPixel{t}{j}, 'dson' );
                    case 4 % ���ѵõ���Сϸ��
                        eval_str = CX_plot_event_sp( eval_str, SuperPixel{t}{j}, 'vson' );
                        
                    case 5 % merge
                        mm = find(Fmj{t}(j,:)==1);
                        source = candidate_k_last{t}{j,mm};
                        SuperPixel{t}{j}.color_index = [ SuperPixel{t-1}{source(1)}.color_index, SuperPixel{t-1}{source(2)}.color_index ];
                        eval_str = CX_plot_event_sp( eval_str, SuperPixel{t}{j}, 'm');   % �ںϵ�����ϸ�����
                        im = CX_fill_color_sp( SuperPixel{t}{j}, im, [1 1 1] );  % �ںϵõ�����Ҫ����ɫ

                    case 6 % appear
                        eval_str = CX_plot_event_sp( eval_str, SuperPixel{t}{j}, 's');  % �³���ϸ�����
                        SuperPixel{t}{j}.color_index = randi(nc); 
                        im = CX_fill_color_sp( SuperPixel{t}{j}, im, color, SuperPixel{t}{j}.color_index ); % ���ֱȽ����⣬�³��ֵ���Ҫ��������ɫ
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
                        SuperPixel{t+1}{e_next}.color_index = SuperPixel{t}{j}.color_index;  

                    case 2 % disappear
                        eval_str = CX_plot_event_sp( eval_str, SuperPixel{t}{j}, 't' );

                    case 3 % divide
                        mm=find(Fid{t}(j,:)==1);
                        % �ҵ�2����ϸ���������Ǽ̳���ɫ�����ñ�ǣ���ɫ�����µģ�
                        son = candidate_k_next{t}{j,mm};

                        SuperPixel{t+1}{son(1)}.color_index = randi(nc);
                        SuperPixel{t+1}{son(2)}.color_index = randi(nc);

                        eval_str = CX_plot_event_sp( eval_str, SuperPixel{t}{j}, 'd' );

                    case 4 % split
                        mm=find(Fiv{t}(j,:)==1);
                        % �ҵ�2������ϸ���������Ǽ̳���ɫ
                        split = candidate_k_next{t}{j,mm};
                        if numel(SuperPixel{t}{j}.color_index)==2 % ����ҵ�2����ɫ�͸�����ȡ������˵����merge��split�ˣ��д���
                            SuperPixel{t+1}{split(1)}.color_index = SuperPixel{t}{j}.color_index(1);
                            SuperPixel{t+1}{split(2)}.color_index = SuperPixel{t}{j}.color_index(2);
                        else
                            SuperPixel{t+1}{split(1)}.color_index = randi(nc);
                            SuperPixel{t+1}{split(2)}.color_index = randi(nc); % �������
                        end
                        eval_str = CX_plot_event_sp( eval_str, SuperPixel{t}{j}, 'v' );  
                        
                    case 5 % merge
                        eval_str = CX_plot_event_sp( eval_str, SuperPixel{t}{j}, 'mson' );
                end
            end
        end
    end

    % ��Ϊim�ڲ��ϸı䣬��˻�ͼ������eval_str�У���im������ٻ�ͼ��2015.12.14���ǣ���
    im = drawregionboundaries(origin, im, [0 0 0]); % ����ԭʼ�ָ���
    imshow(im);hold on;
    eval(eval_str);hold off; 
    
    savename = [output_addr, strrep(png_dir(t).name,'png','fig')];
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

close all


