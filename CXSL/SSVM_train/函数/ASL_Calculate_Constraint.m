function F = ASL_Calculate_Constraint( dataset, s_frame,e_frame,fij,fit,fid,fiv,fmj,fsj )
% ================================================================== %
%
% ���������ڼ���Լ������ F��active sl����ר��
% ������ CXSL_Assign_FlowVar_With_Loss ��Ԥ�������̱����Լ�������ʧ����
% ������������̱�������ʧ�����Լ�Լ������
%
% ASL��Լ����OURS������
%       1��ASL��ÿ��ǰ������һ����˵�������ǵ���˵ǰ�����Ƕ��˵ǰ��
%       2��ASL��ÿ����˵�����뱻����
%
% ================================================================== %
global Ellipse n conflict_fij conflict_pair_next_xy conflict_pair_last_xy;

[ ~, trackpath ] = getpath( dataset );
% �������ݣ�����ѡ������ѵ��������Լ��ϵ�����

if exist('SuperPixel','var') % ����õ��ǳ����أ��Ͳ��ø�ֵ
    disp('  ���õ���SuperPixel��˵��');
    Ellipse = SuperPixel;
    clear SuperPixel
end

F1 = [];
F2 = [];
F3 = [];
F4 = [];

%######################################## ����Լ������ #############################################
disp('   ����Լ��1������Լ������');
assert(s_frame+1==e_frame); % ֻʹ����2֡�����

%% Լ��1
%################################################################
% ������ݼ�3��Fluo-N2DH-SIM+����Լ������ֹmerge��split�¼��ķ�������Բ��˵�У������ز���ֹ��
if ~isempty(strfind(trackpath,'SIM+'))
    disp('Attention! merge&split event has been canceled��')
    for t = s_frame+1:e_frame
        F1 = [ F1, sum(fmj{t}(:))==0 ];
    end
    for t = s_frame:e_frame-1
        F1 = [ F1, sum(fiv{t}(:))==0 ];
    end
end

%% Լ��2
% �ں�֮���������̷���
for t = s_frame+1:e_frame-1
    for j=1:n(t)
        F2 = [ F2, sum(fmj{t}(j,:)) + sum(fid{t}(j,:)) <= 1 ];
    end
end

%% Լ��3��ÿ������˵ǰ�����뱻����
t = e_frame; % 2-t ֡Ҫ��ÿ��ǰ����������һ�������ͣ���ڽ��ͣ�
for j=1:n(t)
    sum_fid = 0;
    sum_fiv = 0;
    sum_fij = 0;
    % ����ǰ�����뱻���ͣ�����ں�Ϊ1
    if Ellipse{t}{j}.num_hypoth == 1
        
        % sum_fij Ϊ����Ǩ�Ƶ� j ��fij֮�ͣ���ںͣ������ںͿ����� sum(fij{t}(j,:)) ��ʾ
        for ind=1:size(conflict_fij{t-1}{j}, 1)
            sum_fij = sum_fij + fij{t-1}( conflict_fij{t-1}{j}(ind,1), conflict_fij{t-1}{j}(ind,2) );
        end
        
        % sum_fid Ϊ���з��ѵ����� j �� pair �� fid ֮��
        for ind=1:numel(conflict_pair_next_xy{t}{j})/2
            sum_fid = sum_fid + fid{t-1}( conflict_pair_next_xy{t}{j}(ind,1), conflict_pair_next_xy{t}{j}(ind,2) );
            sum_fiv = sum_fiv + fiv{t-1}( conflict_pair_next_xy{t}{j}(ind,1), conflict_pair_next_xy{t}{j}(ind,2) );
        end
        F3 = [ F3, sum_fij + fsj{t}(j) + sum(fmj{t}(j,:)) + sum_fid + sum_fiv == 1 ];
    else
        % ��Ŀ��ǰ����������һ��Ҫ������
        %
        % �ⲿ��д��j��j+jplusΪһ��ǰ����
    end
end

% ͬʱ��Ҫ���1~t-1֡�ĳ��ڱ��뱻����
t = s_frame;
for j=1:n(t)
    if Ellipse{t}{j}.num_hypoth == 1
        sum_fmj = 0;
        % sum_fmj Ϊ���а��� j ���ں� pair �� fmj ֮��
        for ind=1:numel(conflict_pair_last_xy{t}{j})/2
            sum_fmj = sum_fmj + fmj{t+1}( conflict_pair_last_xy{t}{j}(ind,1), conflict_pair_last_xy{t}{j}(ind,2) );
        end
        F3 = [ F3, sum(fij{t}(j,:)) + fit{t}(j) + sum(fid{t}(j,:)) + sum(fiv{t}(j,:)) + sum_fmj == 1 ];
    else
        % ��Ŀ��ǰ����������һ��Ҫ������
        %
        % �ⲿ��д��j��j+jplusΪһ��ǰ����
    end
end

disp('   ����Լ��2�����˵Լ������');

%% Լ��4��ÿ�����˵ǰ���н���һ������

for t = s_frame:e_frame
    j = 1;
    while j<=n(t)    % j����
        %% ��ǰ����˵��j����1������Ѱ���¸����˵ǰ��
        if Ellipse{t}{j}.num_hypoth==1
            j = j + 1;
            continue;
        end
        
        % ��һ�������j��j+jplus������ͬһǰ��
        jplus = Ellipse{t}{j}.num_hypoth - 1;
        % �������ʵ����ͬ�Ĺ��ܣ������븴��
        %         jplus = 1;
        %         % j+jplus���ܳ�������
        %         while j+jplus<=numel(Ellipse{t}) && Ellipse{t}{j+jplus}.ind_region == Ellipse{t}{j}.ind_region
        %             jplus = jplus + 1;
        %         end
        %         jplus = jplus - 1;
        %##################### j��j+jplusΪһ��ǰ�� ####################

        %% �˴����һ����Ŀ��ǰ������Ҫ��һ����ڵ�����
        %##########################################
        % ��һ֡�Ķ�Ŀ��ǰ������Ϊ1
        if t==s_frame
            % �ۼ�����
            sum_out = 0;
            for thisj=j:j+jplus
                sum_fmj = 0;
                % ���ǰ��������ϸ��merge���ڵĺ�
                for ind=1:numel(conflict_pair_last_xy{t}{thisj})/2
                    sum_fmj = sum_fmj + fmj{t+1}( conflict_pair_last_xy{t}{thisj}(ind,1), conflict_pair_last_xy{t}{thisj}(ind,2) );
                end
                % ���ǰ��������ϸ������ 4 �ֳ��ڵĺ�
                sum_fij = sum(fij{t}(thisj,:));
                sum_fit = fit{t}(thisj);
                sum_fid = sum(fid{t}(thisj,:));
                sum_fiv = sum(fiv{t}(thisj,:));
                
                sum_out = sum_out + sum_fmj + sum_fij + sum_fit + sum_fid + sum_fiv;
                
            end
            F4 = [ F4,  sum_out==1 ];
        else
            % �ڶ�֡�Ķ�Ŀ��ǰ�����Ϊ1
            sum_in = 0;
            for thisj=j:j+jplus
                
                sum_fij = 0;
                sum_fid = 0;
                sum_fiv = 0;
                % ���ǰ��������ϸ�� split �� divide ��ڵĺ�
                for ind=1:numel(conflict_pair_next_xy{t}{thisj})/2
                    sum_fid = sum_fid + fid{t-1}( conflict_pair_next_xy{t}{thisj}(ind,1), conflict_pair_next_xy{t}{thisj}(ind,2) );
                    sum_fiv = sum_fiv + fiv{t-1}( conflict_pair_next_xy{t}{thisj}(ind,1), conflict_pair_next_xy{t}{thisj}(ind,2) );
                end
                % ���ǰ��������ϸ������ 3 ����ڵĺ�

                % sum_fij Ϊ����Ǩ�Ƶ� j ��fij֮�ͣ���ںͣ������ںͿ����� sum(fij{t}(j,:)) ��ʾ
                for ind=1:size(conflict_fij{t-1}{thisj}, 1)
                    sum_fij = sum_fij + fij{t-1}( conflict_fij{t-1}{thisj}(ind,1), conflict_fij{t-1}{thisj}(ind,2) );
                end

                sum_fmj = sum(fmj{t}(thisj,:));
                sum_fsj = fsj{t}(thisj);
                
                sum_in = sum_in + sum_fmj + sum_fij + sum_fsj + sum_fid + sum_fiv;
                
            end
            F4 = [ F4,  sum_in==1 ];
        end
        %##########################################
        %
        % ����һ������һ����Ŀ��ǰ���������������Ϊ1��Լ��
        % ����֤ÿ����Ŀ��ǰ�������뱻���ͣ��������龰�Ĵ���
        %
        %##############################################################
       j = j + jplus + 1; % jֱ���߳����ǰ�� 
       
    end  %## end while
end

%% �������յ�Ŀ�꺯������Լ��

%#################################
% ����Լ�������Ľ������£�
%######### ��ҪԼ�� 4 �� ##########

% F1��������ݼ�3��Fluo-N2DH-SIM+����Լ������ֹmerge��split�¼��ķ�����
% F2���ں�֮���������̷���
% F3��ÿ������˵ǰ�����뱻����
% F4��ÿ�����˵ǰ���н���һ������


F = [ F1, F2, F3, F4 ];

end



