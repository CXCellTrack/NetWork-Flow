%######################################
%
% 2015.5.30 CX on desk
% ���ã�����������ڽ�������ɵ�ϸ����ɫ�����ڹ۲�
% ���ݴ洢������ɫ��ͼƬ�����ڡ����ӻ����ٱ�ǡ�
% ������ϵ������ CX_RePlot_Ellipse ���Ƹ��º����Բ
%          ���� CX_plot �� CX_plot_event ���и��ٱ��
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

load([trackpath, '\Pair\Pre_data_New.mat']);
fig_addr = [trackpath, '\�����ͼ\'];
 
if 1
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
color_numbel = numel(color)/3;
color = color(randperm(color_numbel),:);

% �����ͼ�װ壬�������ͼ
fig_dir = dir([ fig_addr, '*.fig' ]);
% ���ñ���ͼƬ���ļ��У���������е�figͼƬ
delete([output_addr,'*.fig']);
% ���������ͼ�ļ����ڵ�ͼƬ������Ŀǰ��frameСʱ����Ҫ�ػ����ͼƬ
if 0
    CX_RePlot_Ellipse( dataset );
end

%% ��ʼ����ѭ������
for t = s_frame:e_frame
    fig_name = [ fig_addr, fig_dir(t).name ];
    openfig(fig_name, 'new'); %'invisible'
    disp(['  ���ڴ���', fig_name, '...']);
    hold;
        %##### ��һ֡��Ҫ���������ɫ #####
        if t==s_frame
            for j=1:n(t)
                %#####��һ�δ��뷴���õ���ֱ�Ӵ��#####
                sum_Fmj = 0;
                % sum_fmj Ϊ���а��� j ���ں� pair �� fmj ֮��
                for ind=1:numel(conflict_pair_last_xy{t}{j})/2
                    sum_Fmj = sum_Fmj + Fmj{t+1}( conflict_pair_last_xy{t}{j}(ind,1), conflict_pair_last_xy{t}{j}(ind,2) );
                end
                %####################################
                if sum(Fij{t}(j,:)) + Fit{t}(j) + sum(Fid{t}(j,:)) + sum(Fiv{t}(j,:)) + sum_Fmj ~=0   %%�г��ڵ�ϸ���ͷ�����ɫ
                    Ellipse{t}{j}.color_index = randi(64);
                end
            end
        end
        %###############################
    for j=1:n(t)
        %#################### ��ɫ #######################
        if isfield( Ellipse{t}{j}, 'color_index' )
            % �����Բ��2����ɫ����˵����Ϊmerge��Ǩ�����ģ�������ɫ
            % ע�⣬���ﲻ�ܽ�merge������ɫ����Ϊ���е�����ʱ��û�ж��Ƿ�Ϊmerge
            % merge��ɫ�Ĳ������ж�merge��
            if numel( Ellipse{t}{j}.color_index )==2
                CX_plot( Ellipse{t}{j}, [1 1 1] );
            elseif numel( Ellipse{t}{j}.color_index )==1
                CX_plot( Ellipse{t}{j}, color(Ellipse{t}{j}.color_index,:) );
            end
        end
    
        %################# ���ݲ�ͬ���¼�����ɫ���ݸ���� ##################
        
        %########## ����¼��������ڵ�һ֡ ##########
        if t~=s_frame
            %########## ���� ##########
            if Fsj{t}(j)==1   
                indc = randi(color_numbel);
                CX_plot_event( Ellipse{t}{j}, 's');  % �³���ϸ�����
                Ellipse{t}{j}.color_index = indc; 
                CX_plot( Ellipse{t}{j}, color(Ellipse{t}{j}.color_index,:) ); % ���ֱȽ����⣬�³��ֵ���Ҫ��������ɫ
            end
            %########## �ں� ##########
            if sum(Fmj{t}(j,:))==1
                mm = find(Fmj{t}(j,:)==1);
                source = candidate_k_last{t}{j,mm};
                Ellipse{t}{j}.color_index = [ Ellipse{t-1}{source(1)}.color_index, Ellipse{t-1}{source(2)}.color_index ];
                CX_plot_event( Ellipse{t}{j}, 'm');   % �ںϵ�����ϸ�����
                CX_plot( Ellipse{t}{j}, [1 1 1] );  % �ںϵõ�����Ҫ����ɫ
            end
        end
        
        %########## �����¼������������һ֡ ##########
        if t~=e_frame
            %########### ���� ##########
            if sum(Fid{t}(j,:))==1    
                mm=find(Fid{t}(j,:)==1);
                % �ҵ�2����ϸ���������Ǽ̳���ɫ
                son = candidate_k_next{t}{j,mm};
                
                Ellipse{t+1}{son(1)}.color_index = Ellipse{t}{j}.color_index;
                Ellipse{t+1}{son(2)}.color_index = Ellipse{t}{j}.color_index;

                CX_plot_event( Ellipse{t}{j}, 'd' );
            end

            %########### ���� ##########
            if sum(Fiv{t}(j,:))==1     
                mm=find(Fiv{t}(j,:)==1);
                % �ҵ�2������ϸ���������Ǽ̳���ɫ
                split = candidate_k_next{t}{j,mm};
                if numel(Ellipse{t}{j}.color_index)==2 % ����ҵ�2����ɫ�͸�����ȡ������ͼ̳�
                    Ellipse{t+1}{split(1)}.color_index = Ellipse{t}{j}.color_index(1);
                    Ellipse{t+1}{split(2)}.color_index = Ellipse{t}{j}.color_index(2);
                else
                     Ellipse{t+1}{split(1)}.color_index = Ellipse{t}{j}.color_index;
                     Ellipse{t+1}{split(2)}.color_index = Ellipse{t}{j}.color_index;
                end
                CX_plot_event( Ellipse{t}{j}, 'v' );  
            end
            %########### ��ʧ ##########
            if Fit{t}(j)==1     
                CX_plot_event( Ellipse{t}{j}, 't' );
            end

            %########### Ǩ�� ##########
            if sum(Fij{t}(j,:))==1
                mm = find(Fij{t}(j,:)==1);
                e_next = candidate_fij{t}(j,mm);
                Ellipse{t+1}{e_next}.color_index = Ellipse{t}{j}.color_index;
            end

        end
        
    end
    hold;
    savename = [output_addr, fig_dir(t).name];
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