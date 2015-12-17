function [ Ffull Fbase ] = CXSL_Calculate_Constraint_New_Conflict( dataset, use_op_cons, s_frame,e_frame,fij,fit,fid,fiv,fmj,fsj )
% ================================================================== %
%
% ���������ڼ���Լ������ F
% ������ CXSL_Assign_FlowVar_With_Loss ��Ԥ�������̱����Լ�������ʧ����
% ������������̱�������ʧ�����Լ�Լ������
%
% ================================================================== %

[ ~, trackpath ] = getpath( dataset );
% �������ݣ�����ѡ������ѵ��������Լ��ϵ�����
Pre_data_addr = [ trackpath, '\Pair\Pre_data_New.mat' ];  
load( Pre_data_addr);
if exist('SuperPixel','var') % ����õ��ǳ����أ��Ͳ��ø�ֵ
    disp('  ���õ���SuperPixel��˵��');
    Ellipse = SuperPixel;
    clear SuperPixel
end
    
F1 = [];
F2 = [];
Fop1 = [];
Fop2 = [];
Fop3 = [];
Fop4 = [];
Fop5 = [];
F3 = [];
F4 = [];
Fop6 = [];
Fop7 = [];

%######################################## ����Լ������ #############################################
disp('   ����Լ��1�������غ㡭��');
tic
%% Լ������1�� �����غ�

if e_frame - s_frame ~= 1 % ������֡����Ϊ2ʱ���д�����
    disp('      1)һ����Լ������');
    %#############  �м�֡����ڳ����غ�  ###########
    for t = s_frame+1:e_frame-1
        disp(['      ���ڼ����',num2str(t),'֡�Ľ����غ�...']);
        for j=1:n(t)
            sum_fij = 0;
            sum_fid = 0;
            sum_fiv = 0;
            sum_fmj = 0;

            % sum_fij Ϊ����Ǩ�Ƶ� j ��fij֮�ͣ���ںͣ������ںͿ����� sum(fij{t}(j,:)) ��ʾ 
            for ind=1:size(conflict_fij{t-1}{j}, 1)
                sum_fij = sum_fij + fij{t-1}( conflict_fij{t-1}{j}(ind,1), conflict_fij{t-1}{j}(ind,2) );
            end

            % sum_fid Ϊ���з��ѵ����� j �� pair �� fid ֮��
            for ind=1:numel(conflict_pair_next_xy{t}{j})/2
                sum_fid = sum_fid + fid{t-1}( conflict_pair_next_xy{t}{j}(ind,1), conflict_pair_next_xy{t}{j}(ind,2) );
                sum_fiv = sum_fiv + fiv{t-1}( conflict_pair_next_xy{t}{j}(ind,1), conflict_pair_next_xy{t}{j}(ind,2) );
            end
            % sum_fmj Ϊ���а��� j ���ں� pair �� fmj ֮��
            for ind=1:numel(conflict_pair_last_xy{t}{j})/2
                sum_fmj = sum_fmj + fmj{t+1}( conflict_pair_last_xy{t}{j}(ind,1), conflict_pair_last_xy{t}{j}(ind,2) );
            end
            F1 = [ F1, sum_fij + fsj{t}(j) + sum(fmj{t}(j,:)) + sum_fid + sum_fiv == sum(fij{t}(j,:)) + fit{t}(j) + sum(fid{t}(j,:)) + sum(fiv{t}(j,:)) + sum_fmj ];
        end
    end
    
end
% end if

%#############  2��t֡���Ψһ  ###########
disp('      2)���Լ������');
for t = s_frame+1:e_frame
    for j=1:n(t)
        sum_fid = 0;
        sum_fiv = 0;
        sum_fij = 0;
        % sum_fij Ϊ����Ǩ�Ƶ� j ��fij֮�ͣ���ںͣ������ںͿ����� sum(fij{t}(j,:)) ��ʾ 
        for ind=1:size(conflict_fij{t-1}{j}, 1)
            sum_fij = sum_fij + fij{t-1}( conflict_fij{t-1}{j}(ind,1), conflict_fij{t-1}{j}(ind,2) );
        end
        
        % sum_fid Ϊ���з��ѵ����� j �� pair �� fid ֮��
        for ind=1:numel(conflict_pair_next_xy{t}{j})/2
            sum_fid = sum_fid + fid{t-1}( conflict_pair_next_xy{t}{j}(ind,1), conflict_pair_next_xy{t}{j}(ind,2) );
            sum_fiv = sum_fiv + fiv{t-1}( conflict_pair_next_xy{t}{j}(ind,1), conflict_pair_next_xy{t}{j}(ind,2) );
        end
        % ��������һ��   Ǩ����           �³���      merge�ں���      ��ϸ��
        F2 = [ F2, sum_fij + fsj{t}(j) + sum(fmj{t}(j,:)) + sum_fid + sum_fiv <= 1 ];
    end
end

%### 1��t-1֡����Ψһ����һ����Լ������ʱ���ɼ�Ϊ����һ֡����Ψһ��  ###%
disp('      3)����Լ������');
for t=s_frame%:e_frame-1
    for j=1:n(t)
        sum_fmj = 0;
        % sum_fmj Ϊ���а��� j ���ں� pair �� fmj ֮��
        for ind=1:numel(conflict_pair_last_xy{t}{j})/2
            sum_fmj = sum_fmj + fmj{t+1}( conflict_pair_last_xy{t}{j}(ind,1), conflict_pair_last_xy{t}{j}(ind,2) );
        end
        % ���������һ��   Ǩ�Ƴ�ȥ        ��ʧ        ĸϸ��        ��һ֡split        ��һ֡�ں�
        F2 = [ F2, sum(fij{t}(j,:)) + fit{t}(j) + sum(fid{t}(j,:)) + sum(fiv{t}(j,:)) + sum_fmj <= 1 ];
    end
end

toc;
disp('   ����Լ��2����ѡԼ������');
disp(use_op_cons);
tic

%% ��ѡԼ��1
%################################################################
% ������ݼ�3��Fluo-N2DH-SIM+����Լ������ֹmerge��split�¼��ķ�������Բ��˵�У������ز���ֹ��
% if ~isempty(strfind(trackpath,'SIM+'))
%     disp('Attention! merge&split event has been canceled��')
%     for t = s_frame+1:e_frame
%         Fop1 = [ Fop1, sum(fmj{t}(:))<=0 ];
%     end
%     for t = s_frame:e_frame-1
%         Fop1 = [ Fop1, sum(fiv{t}(:))<=0 ];
%     end
% end

if any(use_op_cons==1)
    %################## ��ѡ���Լ��һ #################������ʵ����е�ì�ܣ�
    % ���ѳ�ȥ�Ĳ�������ͬһ��ǰ����
    for t = s_frame:e_frame-1
        for j=1:n(t)
            tmp_conflict = 0;
            for mm=1:6
                sons = candidate_k_next{t}{j,mm};
                if Ellipse{t+1}{sons(1)}.ind_region == Ellipse{t+1}{sons(2)}.ind_region % ���ѳ�ȥ�Ĳ�������ͬһ��ǰ����
                    tmp_conflict = tmp_conflict + fid{t}(j,mm);
                end
            end
            % ��������ڶ��ǵ���ǰ����tmp_conflict����Ϊdouble 0�����Ҫ�ж�
            if ~isa(tmp_conflict, 'double')
                Fop1 = [ Fop1, tmp_conflict <= 0];
            end
        end
    end
end

%% ��ѡԼ��2
if any(use_op_cons==2)
    %################## ��ѡ���Լ���� #################������ʵ����е�ì�ܣ�
    % ��������ͬһ��ǰ����ϸ�������ں� ��pair��ʽ��ʵ����������ף�ֻ��Ҫ�Ƚ��ں�pair�е�2ϸ���Ƿ���һ��ǰ���ھͿ�����
    for t = s_frame+1:e_frame
        for j=1:n(t)
            tmp_conflict = 0;
            for mm=1:6
                sources = candidate_k_last{t}{j,mm};
                if Ellipse{t-1}{sources(1)}.ind_region ~= Ellipse{t-1}{sources(2)}.ind_region 
                    tmp_conflict = tmp_conflict + fmj{t}(j,mm);
                end
            end
            % ��������ڶ��ǵ���ǰ����tmp_conflict����Ϊ0�����Ҫ�ж�
            if ~isa(tmp_conflict, 'double')
                Fop2 = [ Fop2, tmp_conflict <= 0];
            end
        end
    end
end

%% ��ѡԼ��3
if any(use_op_cons==3)
    % %################## ��ѡ���Լ���� #################
    % �ں�֮���������̷���
    for t = s_frame+1:e_frame-1
        for j=1:n(t)
            Fop3 = [ Fop3, sum(fmj{t}(j,:)) + sum(fid{t}(j,:)) <= 1 ];
        end
    end
end

%% ��ѡԼ��4�Ǵ���ģ��������ã�
    %################## ��ѡ���Լ���� #################�������Լ����ע�͵���
    % ��Ϊ����˵ǰ�����Ǳ��뱻���ͣ������˵ǰ����Լ��д��ì�ܼ��ϴ�����˴˴�ע�͵�
    %     for t = s_frame+1:e_frame % 2-t ֡Ҫ��ÿ��ǰ����������һ�������ͣ���ڽ��ͣ�
    %         for j=1:n(t)
    %             sum_fid = 0;
    %             sum_fiv = 0;
    %             sum_fij = 0;
    %             % ����ǰ�����뱻���ͣ�����ں�Ϊ1
    %             if Ellipse{t}{j}.num_hypoth == 1
    % 
    %                 % sum_fij Ϊ����Ǩ�Ƶ� j ��fij֮�ͣ���ںͣ������ںͿ����� sum(fij{t}(j,:)) ��ʾ 
    %                 for ind=1:size(conflict_fij{t-1}{j}, 1)
    %                     sum_fij = sum_fij + fij{t-1}( conflict_fij{t-1}{j}(ind,1), conflict_fij{t-1}{j}(ind,2) );
    %                 end
    % 
    %                 % sum_fid Ϊ���з��ѵ����� j �� pair �� fid ֮��
    %                 for ind=1:numel(conflict_pair_next_xy{t}{j})/2
    %                     sum_fid = sum_fid + fid{t-1}( conflict_pair_next_xy{t}{j}(ind,1), conflict_pair_next_xy{t}{j}(ind,2) );
    %                     sum_fiv = sum_fiv + fiv{t-1}( conflict_pair_next_xy{t}{j}(ind,1), conflict_pair_next_xy{t}{j}(ind,2) );
    %                 end
    %                 Fop4 = [ Fop4, sum_fij + fsj{t}(j) + sum(fmj{t}(j,:)) + sum_fid + sum_fiv >= 1 ];
    %             else
    %                 % ��Ŀ��ǰ����������һ��Ҫ������
    %                 %
    %                 % �ⲿ��д��j��j+jplusΪһ��ǰ���� 
    %             end
    %         end
    %     end
    % 
    %     t = s_frame; % ͬʱ��Ҫ��ӵ�һ֡�ĳ��ڱ��뱻����
    %     for j=1:n(t)
    %         if Ellipse{t}{j}.num_hypoth == 1
    %             sum_fmj = 0;
    %             % sum_fmj Ϊ���а��� j ���ں� pair �� fmj ֮��
    %             for ind=1:numel(conflict_pair_last_xy{t}{j})/2
    %                 sum_fmj = sum_fmj + fmj{t+1}( conflict_pair_last_xy{t}{j}(ind,1), conflict_pair_last_xy{t}{j}(ind,2) );
    %             end
    %             Fop4 = [ Fop4, sum(fij{t}(j,:)) + fit{t}(j) + sum(fid{t}(j,:)) + sum(fiv{t}(j,:)) + sum_fmj >= 1 ];
    %         else
    %             % ��Ŀ��ǰ����������һ��Ҫ������
    %             %
    %             % �ⲿ��д��j��j+jplusΪһ��ǰ����
    %         end
    %     end

%% ��ѡ��Լ��5 
if any(use_op_cons==5)
    %##################  #################
    % Ҫ�����split��ȥ��2����ϸ��������ͬһǰ����
    for t = s_frame:e_frame-1
        for j=1:n(t)
            tmp_conflict = 0;
            for mm=1:6
                sons = candidate_k_next{t}{j,mm};
                if Ellipse{t+1}{sons(1)}.ind_region ~= Ellipse{t+1}{sons(2)}.ind_region % ���ѳ�ȥ�ı����ǵ���ǰ��
                    tmp_conflict = tmp_conflict + fiv{t}(j,mm);
                end
            end
            % ��������ڶ��ǵ���ǰ����tmp_conflict����Ϊdouble 0�����Ҫ�ж�
            if ~isa(tmp_conflict, 'double')
                Fop5 = [ Fop5, tmp_conflict <= 0];
            end
        end
    end
    %#################################################
end
toc
disp('   ����Լ��3��ì��Լ������');
tic

%% Լ������3�� ì�ܼ�˵�����ų� 

conflict = cell(e_frame, 1);
for t = s_frame:e_frame-1
    conflict{t}={};
    j = 1;
    while j<=n(t)    % j����
        % ��ǰ����˵��j����1������Ѱ���¸����˵ǰ��
        if Ellipse{t}{j}.num_hypoth==1
            j = j + 1;
            continue;
        end
        
        n_basic = numel(Ellipse{t}{j}.flag_combine);
        if n_basic==1
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
            
        %% A. �˴����һ����Ŀ��ǰ������Ҫ��һ����ڵ�����
        if 0
            %##########################################
            % ��һ֡�Ķ�Ŀ��ǰ����������Ϊ1���Ǳ�ҪԼ����
            if t==s_frame
                sum_fmj = 0;
                sum_fij = 0;
                sum_fit = 0;
                sum_fid = 0;
                sum_fiv = 0;
                for j_muti=j:j+jplus
                    % ���ǰ��������ϸ��merge���ڵĺ�
                    for ind=1:numel(conflict_pair_last_xy{t}{j_muti})/2
                        sum_fmj = sum_fmj + fmj{t+1}( conflict_pair_last_xy{t}{j_muti}(ind,1), conflict_pair_last_xy{t}{j_muti}(ind,2) );
                    end
                    % ���ǰ��������ϸ������ 4 �ֳ��ڵĺ�
                    sum_fij = sum_fij + sum(fij{t}(j_muti,:));
                    sum_fit = sum_fit + fit{t}(j_muti);
                    sum_fid = sum_fid + sum(fid{t}(j_muti,:));
                    sum_fiv = sum_fiv + sum(fiv{t}(j_muti,:));
                end
                Fop6 = [ Fop6, sum_fmj + sum_fij + sum_fit + sum_fid + sum_fiv >= 1 ];

            else
            % 2-t֡�Ķ�Ŀ��ǰ���������Ϊ1���Ǳ�ҪԼ����
                sum_fmj = 0;
                sum_fij = 0;
                sum_fsj = 0;
                sum_fid = 0;
                sum_fiv = 0;
                for j_muti=j:j+jplus
                    % ���ǰ��������ϸ�� split �� divide ��ڵĺ�
                    for ind=1:numel(conflict_pair_next_xy{t}{j_muti})/2
                        sum_fid = sum_fid + fid{t-1}( conflict_pair_next_xy{t}{j_muti}(ind,1), conflict_pair_next_xy{t}{j_muti}(ind,2) );
                        sum_fiv = sum_fiv + fiv{t-1}( conflict_pair_next_xy{t}{j_muti}(ind,1), conflict_pair_next_xy{t}{j_muti}(ind,2) );
                    end
                    % ���ǰ��������ϸ������ 3 ����ڵĺ�

                    sum_fij_m = 0;
                    % sum_fij Ϊ����Ǩ�Ƶ� j ��fij֮�ͣ���ںͣ������ںͿ����� sum(fij{t}(j,:)) ��ʾ 
                    for ind=1:size(conflict_fij{t-1}{j_muti}, 1)
                        sum_fij_m = sum_fij_m + fij{t-1}( conflict_fij{t-1}{j_muti}(ind,1), conflict_fij{t-1}{j_muti}(ind,2) );
                    end

                    % �� sum_fij_m ��ںʹ���ע�͵��� sum(fij{t-1}(:,j_muti))
                    sum_fij = sum_fij + sum_fij_m;
    %                     sum_fij = sum_fij + sum(fij{t-1}(:,j_muti));
                    sum_fmj = sum_fmj + sum(fmj{t}(j_muti,:));
                    sum_fsj = sum_fsj + fsj{t}(j_muti);
                end
                Fop7 = [ Fop7, sum_fmj + sum_fij + sum_fsj + sum_fid + sum_fiv >= 1 ];
            end
        end
            %##########################################   
            %
            % ����һ������һ����Ŀ��ǰ���������������Ϊ1��Լ��
            % ����֤ÿ����Ŀ��ǰ�������뱻���ͣ��������龰�Ĵ���
            %
            %##############################################################

        %% B. conflict��ʾ���п��ܵ�ì�ܼ���
        use_22_cons = 0;
        if use_22_cons % �·�����22ì��Լ��
            % ----------------------------------------------------------- %
            % �·�������ì��Լ��д��һ���������ҳ�ǰ���ڵ�22ì�ܣ���һ��������
            conflict_New = [];
            cc = 1; % ���ڼ���
            for uu=j:j+jplus-1
                for vv=uu+1:j+jplus
                    if ~isequal( Ellipse{t}{uu}.flag_combine & Ellipse{t}{vv}.flag_combine, zeros(1,n_basic)) % ���벻Ϊ0�������ì�ܼ�
                        conflict_New(cc,:) = [uu vv]; % ÿһ����һ������ì��
                        cc = cc + 1;
                    end
                end
            end
            % �����ⲿ�����������conflict�е�ì��Լ����conflict���������ڸ�Լ���ĸ����������Դ���֤Լ������ȷ�ԣ�
            [ F_j_jlus ] = Calculate_Conflict_Constraint( conflict_New, t, s_frame, fij, fit, fid, fiv, fmj, fsj,...
                conflict_pair_last_xy, conflict_fij, conflict_pair_next_xy);

            F3 = [ F3, F_j_jlus ];

        end
        % ----------------------------------------------------------- %
        if ~use_22_cons % �ɷ�����Щ���󣬲����ܰ������е�ì�ܼ��ϣ���һ�е�for�������⣩ 2015.7.7���޸���ͨ���� CX_Ellipse_Optimal �м�������

            danduan_flag = eye(n_basic); % ��λ��

            for uu=1:n_basic % ������λ���ÿһ��
                conflict{t}{j,uu} = [];
                for vv=j:j+jplus % ����һ��ǰ��
                    % ���벻Ϊ0�������ì�ܼ�
                    if ~isequal( danduan_flag(uu,:) & Ellipse{t}{vv}.flag_combine, zeros(1,n_basic)) 
                        conflict{t}{j,uu} = [conflict{t}{j,uu}, vv];
                    end
                end

                if numel(conflict{t}{j,uu})==1 % ���é�ܼ���ֻ��1����Բ��������
                    continue;
                end 
                %##########################�˴���������Լ��###############
                conflict_set = conflict{t}{j,uu}; % ì�ܼ��ϵ���ʱ����
                if t==s_frame % ��һ֡���ǳ���
                    % ----------------------------------------------- %
                    conflict_sum_out = 0;
                    %### ��һ֡һ��ì�ܼ��ڳ�������Ϊ1��
                    for jj=1:numel(conflict_set) % conflict_tmp_setΪ[21, 23]����ʽ
                        this_j = conflict_set(jj); % ȡ��ì�ܼ���һ��ϸ���ı��
                        sum_fmj = 0;
                        % sum_fmj Ϊ�ںϳ��ں�
                        tmp_xy = conflict_pair_last_xy{t}{this_j}; % ȡ��������ǰϸ�� ind_j ��pair�������
                        for ind=1:numel(tmp_xy)/2
                            sum_fmj = sum_fmj + fmj{t+1}( tmp_xy(ind,1), tmp_xy(ind,2) );
                        end
                        % �������ں�
                        this_out = sum(fij{t}(this_j,:))+fit{t}(this_j)+sum(fid{t}(this_j,:)) + sum(fiv{t}(this_j,:))+sum_fmj;
                        % ì�ܼ�ȫ����ں�
                        conflict_sum_out = conflict_sum_out + this_out;
                    end
                    F3 = [ F3, conflict_sum_out <= 1 ];   % ����Ψһ
                    % ----------------------------------------------- %
                else
                    % ----------------------------------------------- %
                    % 2��t֡�������
                    conflict_sum_in = 0;
                    %### һ��ì�ܼ����������Ϊ1��
                    for jj=1:numel(conflict_set)
                        this_j = conflict_set(jj); % ȡ��ì�ܼ��е�һ��ϸ���ı��
                        
                        sum_fid = 0;
                        sum_fiv = 0;
                        sum_fij = 0;
                        % =========================================== %
                        % sum_fij Ϊ����Ǩ�Ƶ� j ��fij֮�ͣ���ںͣ������ںͿ����� sum(fij{t}(j,:)) ��ʾ 
                        tmp_xy = conflict_fij{t-1}{this_j};
                        for ind=1:size(tmp_xy, 1)
                            sum_fij = sum_fij + fij{t-1}(tmp_xy(ind,1), tmp_xy(ind,2));
                        end
                        % =========================================== %
                        % sum_fid Ϊ���з��ѳ��ں�
                        tmp_xy = conflict_pair_next_xy{t}{this_j}; % ȡ��������ǰϸ�� ind_j ��pair�������
                        for ind=1:size(tmp_xy,1)
                            sum_fid = sum_fid + fid{t-1}( tmp_xy(ind,1), tmp_xy(ind,2) );
                            sum_fiv = sum_fiv + fiv{t-1}( tmp_xy(ind,1), tmp_xy(ind,2) );
                        end
                        % ������ں�
                        this_in = fsj{t}(this_j)+sum(fmj{t}(this_j,:))+sum_fij+sum_fid+sum_fiv;
                        % ì�ܼ�ȫ����ں�
                        conflict_sum_in = conflict_sum_in + this_in;
                    end    
                    F4 = [ F4, conflict_sum_in <= 1 ]; %%2��n֡����ں�Ψһ
                    % ----------------------------------------------- %
                end
                %##########################################################
            end

        end % end if 0
        j = j + jplus + 1; % jֱ���߳����ǰ��
            
    end  %## end while
end

toc

%% �������յ�Ŀ�꺯������Լ��

%#################################
% ����Լ�������Ľ������£�
%######### ��ҪԼ�� 4 �� ##########

% F1��2��t-1 ֡�� ÿ����Բ�� input = output
% F2��2��t   ֡�� ÿ����Բ�� input<=1 && �� 1 ֡�� output <=1
% F3��1      ֡�� ÿ��ì�ܼ��� output<=1
% F4��2��t   ֡�� ÿ��ì�ܼ��� inputput<=1

%######### ��ѡԼ�� 6 �� ##########

% Fop1��2��t   ֡�� ���ѵõ�����ϸ��������ͬһ��ǰ���� 
% Fop2��1��t-1 ֡�� ����ͬһ��ǰ����ϸ�����ܷ����ں� 
% Fop3��2��t-1 ֡�� �ںϵ����Ĵ�ϸ�����������̷���
% Fop4��1��t   ֡�� ����˵ǰ�����뱻���ͣ����� 2��t ֡����ڽ��ͺ͵� 1 ֡�ĳ��ڽ��ͣ������Լ���Ǵ���ģ���
% Fop5��2��t-1 ֡�� ����õ���2��ϸ��������ͬһ��ǰ���� 
% Fop6��1      ֡�� ���˵ǰ�� ��������Ϊ1
% Fop7��2��t   ֡�� ���˵ǰ�� �������Ϊ1

%##################################
%##########################################################################
% sdpsettings:      �����ѡ����������
% http://users.isy.liu.se/johanl/yalmip/pmwiki.php?n=Commands.Sdpsettings

% verbose           ��ӡ��Ϣ��
% solver            ָ�������
% saveduals         ����Ϊ0���������ż�����Խ�ʡ�ڴ�
% usex0             ��� assign �������������ó�ʼֵ���Լӿ�����ٶ�
% assign            http://users.isy.liu.se/johanl/yalmip/pmwiki.php?n=Commands.Assign
% savesdpafile      ����Լ��������Ŀ�꺯�� http://users.isy.liu.se/johanl/yalmip/pmwiki.php?n=Commands.Savesdpafile
% loadsedumidata    ����Լ��������Ŀ�꺯�� http://users.isy.liu.se/johanl/yalmip/pmwiki.php?n=Commands.Loadsedumidata
      
%##########################################################################

Fbase = [ F1, F2, F3, F4 ];
Foptional = [ Fop1 Fop2 Fop3 Fop5 ]; % 4����6��7�������ã�1��2����ʵ�����ì�ܣ����ֻ����3��5��2015.10.3��
% clear Fop1 Fop2 Fop3 Fop4 Fop5 Fop6 Fop7;
Ffull = [ Fbase, Foptional ];


% clear F1 F2 F3 F4 Foptional;
% ����Ŀ�꺯����Լ����������ջ���
% clearvars -except frame Ellipse F object_function fij fid fiv fit fsj fmj s_frame e_frame fai_x_z sum_cost; 

end



function [ F_j_jlus ] = Calculate_Conflict_Constraint( conflict, t, s_frame, fij, fit, fid, fiv, fmj, fsj,...
    conflict_pair_last_xy, conflict_fij, conflict_pair_next_xy)
% 
% 

F_j_jlus = []; % ������ʾj��jplus��һ��ǰ���ڵ�ì��

%##########################�˴���������Լ��###############
if t==s_frame % ��һ֡���ǳ���
    
    %### ��һ֡2��ì��ϸ���ĳ��ں�����Ϊ1
    for i=1:size(conflict,1)
        sum_fmj = 0;
        j = conflict(i,1); % ȡ��ì�ܼ��е�һ��ϸ���ı��
        k = conflict(i,2);
        constrain_sum1 = sum( fij{t}(k,:) ) + sum( fij{t}(j,:) ); % Ǩ�Ƴ��ں�
        constrain_sum2 = fit{t}(k) + fit{t}(j);   % ��ʧ���ں�
        constrain_sum3 = sum( fid{t}(k,:) ) + sum( fiv{t}(k,:) ) + sum( fid{t}(j,:) ) + sum( fiv{t}(j,:) );% ����/���ѳ��ں�
        % sum_fmj Ϊ�ںϳ��ں�
        tmp_xy_fmj = unique([ conflict_pair_last_xy{t}{j};conflict_pair_last_xy{t}{k} ], 'rows'); % �ҳ�����ì�����꣬��ȥ���ظ��� 
        for ind=1:size(tmp_xy_fmj,1)
            sum_fmj = sum_fmj + fmj{t+1}( tmp_xy_fmj(ind,1), tmp_xy_fmj(ind,2) );
        end
        F_j_jlus = [ F_j_jlus, constrain_sum1 + constrain_sum2 + constrain_sum3 + sum_fmj <= 1 ];   %%����Ψһ
    end
    % ����return������Ҫ�������
    return;
end

% 2��t֡�������
% 2��ì��ϸ������ں�����Ϊ1
for i=1:size(conflict,1)

    j = conflict(i,1); % ȡ��ì�ܼ��е�һ��ϸ���ı��
    k = conflict(i,2);
    % ----------------------------------------------- %
    sum_fij = 0;
    % sum_fij Ϊ����Ǩ�Ƶ� j ��fij֮�ͣ���ںͣ������ںͿ����� sum(fij{t}(j,:)) ��ʾ 
    tmp_xy_fij =  unique([conflict_fij{t-1}{j};conflict_fij{t-1}{k}], 'rows');
    for ind=1:size(tmp_xy_fij, 1)
        sum_fij = sum_fij + fij{t-1}( tmp_xy_fij(ind,1), tmp_xy_fij(ind,2) );
    end
    % ----------------------------------------------- %
    % �� sum_fij_m ���� sum( fij{t-1}(:,ind_j) )
    constrain_sum4 = sum_fij + fsj{t}(j) + sum(fmj{t}(j,:)) + fsj{t}(k) + sum(fmj{t}(k,:)); % Ǩ�ơ����֡��ں�
    % ----------------------------------------------- %
    % sum_fid Ϊ���з��ѳ��ں�                
    sum_fid = 0;
    sum_fiv = 0;
    tmp_xy_fidv = unique([conflict_pair_next_xy{t}{j};conflict_pair_next_xy{t}{k}], 'rows'); % ȡ��������ǰϸ�� ind_j ��pair�������
    for ind=1:size(tmp_xy_fidv,1)
        sum_fid = sum_fid + fid{t-1}( tmp_xy_fidv(ind,1), tmp_xy_fidv(ind,2) );
        sum_fiv = sum_fiv + fiv{t-1}( tmp_xy_fidv(ind,1), tmp_xy_fidv(ind,2) );
    end
    % ----------------------------------------------- %
    F_j_jlus = [ F_j_jlus, constrain_sum4 + sum_fid + sum_fiv <= 1 ]; %%2��n֡����ں�Ψһ
end    

end






