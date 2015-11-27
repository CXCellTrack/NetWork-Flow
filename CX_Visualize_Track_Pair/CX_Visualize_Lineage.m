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
trackpath = 'E:\datasets\first_edition\training_datasets\N2DL-HeLa\02_4-16_track';
load([trackpath, '\Pair\Pre_data_New.mat']);

disp('  ������ʵ���̱�������...');
load([trackpath, '\GT\GT_after_hand_tune\GT_Flow_Variables_New.mat']);
output_addr = [trackpath, '\GT\'];
s_frame = 1;
e_frame = numel(Fmj);
frame = e_frame - s_frame + 1;

color = colormap(hsv);close('1'); % ѡ����ɫ
nc = numel(color)/3;
color = color(randperm(nc),:);


%% ��ʼ����ѭ������
im_addr = 'E:\datasets\first_edition\training_datasets\N2DL-HeLa\02\t00.tif'; % ��ͼ�װ�
im = imread(im_addr);im = imadjust(im);
im = zeros(size(im));
imshow(im);%imcontrast(gca);
hold on;

for t = s_frame:e_frame
    disp(['  ���ڴ���� ', num2str(t), ' ֡...']);
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
                % ����԰
                CX_plot_ball( Ellipse{t}{j}, color );
            end
        end
    end
        %###############################
    for j=1:n(t)
        %################# ���ݲ�ͬ���¼�����ɫ���ݸ���� ##################
        %########## ����¼��������ڵ�һ֡ ##########
        if t~=s_frame
            %########## ���� ##########
            if Fsj{t}(j)==1   
                indc = randi(nc);
                Ellipse{t}{j}.color_index = indc; 
                CX_plot_ball( Ellipse{t}{j}, color ); % ���ֱȽ����⣬�³��ֵ���Ҫ��������ɫ
            end
            %########## �ں� ##########
            if sum(Fmj{t}(j,:))==1
                mm = find(Fmj{t}(j,:)==1);
                source = candidate_k_last{t}{j,mm};
                
                % ����һ���쳣����� 
                Ellipse{t-1}{source(1)} = add_color_index( Ellipse{t-1}{source(1)} );
                Ellipse{t-1}{source(2)} = add_color_index( Ellipse{t-1}{source(2)} );
                
                Ellipse{t}{j}.color_index = [ Ellipse{t-1}{source(1)}.color_index, Ellipse{t-1}{source(2)}.color_index ];
                % ����Ǩ����
                CX_plot_line( Ellipse{t-1}{source(1)}, Ellipse{t}{j}, color );
                CX_plot_line( Ellipse{t-1}{source(2)}, Ellipse{t}{j}, color );
                CX_plot_ball( Ellipse{t}{j}, [1 1 1] );  % �ںϵõ�����Ҫ�»���
            end
        end
        
        %########## �����¼������������һ֡ ##########
        if t~=e_frame
            %########### ���� ##########
            if sum(Fid{t}(j,:))==1    
                mm=find(Fid{t}(j,:)==1);
                % �ҵ�2����ϸ���������Ǽ̳���ɫ
                son = candidate_k_next{t}{j,mm};
                
                Ellipse{t}{j} = add_color_index(Ellipse{t}{j});
                % ����һ���쳣�����
                if numel(Ellipse{t}{j}.color_index)==2 
                    Ellipse{t+1}{son(1)}.color_index = Ellipse{t}{j}.color_index(1);
                    Ellipse{t+1}{son(2)}.color_index = Ellipse{t}{j}.color_index(2);
                elseif numel(Ellipse{t}{j}.color_index)==1
                    Ellipse{t+1}{son(1)}.color_index = Ellipse{t}{j}.color_index;
                    Ellipse{t+1}{son(2)}.color_index = Ellipse{t}{j}.color_index;
                end

                % ���Ʒ����ߺ�����
                CX_plot_line( Ellipse{t}{j}, Ellipse{t+1}{son(1)}, color );
                CX_plot_ball( Ellipse{t+1}{son(1)}, color );
                CX_plot_line( Ellipse{t}{j}, Ellipse{t+1}{son(2)}, color );
                CX_plot_ball( Ellipse{t+1}{son(2)}, color );
            end

            %########### ���� ##########
            if sum(Fiv{t}(j,:))==1     
                mm=find(Fiv{t}(j,:)==1);
                % �ҵ�2������ϸ���������Ǽ̳���ɫ
                split = candidate_k_next{t}{j,mm};
                
                % ����һ���쳣�����
                Ellipse{t}{j} = add_color_index(Ellipse{t}{j});
                
                if numel(Ellipse{t}{j}.color_index)==2 % ����ҵ�2����ɫ�͸�����ȡ������ͼ̳�
                    Ellipse{t+1}{split(1)}.color_index = Ellipse{t}{j}.color_index(1);
                    Ellipse{t+1}{split(2)}.color_index = Ellipse{t}{j}.color_index(2);
                elseif numel(Ellipse{t}{j}.color_index)==1
                     Ellipse{t+1}{split(1)}.color_index = Ellipse{t}{j}.color_index;
                     Ellipse{t+1}{split(2)}.color_index = Ellipse{t}{j}.color_index;
                end
                % ����split�ߺ�����
                CX_plot_line( Ellipse{t}{j}, Ellipse{t+1}{split(1)}, color );
                CX_plot_ball( Ellipse{t+1}{split(1)}, color );
                CX_plot_line( Ellipse{t}{j}, Ellipse{t+1}{split(2)}, color );
                CX_plot_ball( Ellipse{t+1}{split(2)}, color );  
            end
            %########### ��ʧ ##########
            if Fit{t}(j)==1     
%                 CX_plot_event( Ellipse{t}{j}, 't' );
            end

            %########### Ǩ�� ##########
            if sum(Fij{t}(j,:))==1
                mm = find(Fij{t}(j,:)==1);
                e_next = candidate_fij{t}(j,mm);
                
                % ����һ���쳣�����
                Ellipse{t}{j} = add_color_index(Ellipse{t}{j});
                Ellipse{t+1}{e_next}.color_index = Ellipse{t}{j}.color_index;
                CX_plot_line( Ellipse{t}{j}, Ellipse{t+1}{e_next}, color );
            end

        end
        
    end

end
hold off;
