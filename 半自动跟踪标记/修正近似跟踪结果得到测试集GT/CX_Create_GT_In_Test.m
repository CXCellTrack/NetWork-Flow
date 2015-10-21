% ================================================================== %
%
% CX 2015.7.21
% ����ű����ڸ��� CX_Hand_Annotation_In_Test �õ��Ľ��������ɾ������ӣ�
% ��һ����Ϊ�ӽ� gt �� trackdata ����Ϊ ground truth
% ��Ҫ���������1 ɾȥ����켣
%              2 �����µ���ȷ�켣
%
% ================================================================== %

clear;close all;

%% ������Լ����ϵĸ��ٽ��
if 1
    dataset = 'competition';
else
    dataset = 'training';
end

[ ~, trackpath ] = getpath( dataset );

load([ trackpath, '\Pair\Pre_data_New.mat']);
track_data_addr = [ trackpath, '\�ṹ��ѧϰ\Tracking_Data.mat'];
if strcmp(dataset,'training') % �޸�ѵ�����𰸵�ַ������
    track_data_addr = [ trackpath, '\GT\GT_Flow_Variables_New.mat'];
end

load( track_data_addr );
frame = numel(Fmj);
% frame = 20;
% ��������Ľ��
Fij_c = Fij;
Fit_c = Fit;
Fid_c = Fid;
Fiv_c = Fiv;
Fmj_c = Fmj;
Fsj_c = Fsj;

%% �򿪻�����ɵĸ���ͼƬ������GT���� CX_Hand_Annotation_In_Test �н��в�����
load([ trackpath, '\GT\Hand_GT_New.mat']); 
GT_delete = GT_delete_s;
GT_move = GT_move_s;

%% 1�����ӻ�����µĹ켣������Ҫɾ��1��ȥ·��2����·�������1->2�Ĺ켣
for t=1:numel(GT_move)
    % Ϊ�����������ӿ��ٶ�
    if isempty(GT_move{t})
        continue;
    end
    
    for ind=1:size(GT_move{t},1)
        % �ҵ�һ�м�¼������
        rowdata = GT_move{t}(ind,:);
        if isequal(rowdata, [0 0 0 0]) % ��һ��Ϊ�ռ�¼������
            break;
        end

        %% A ɾ��Ҫ�޸ĵĹ켣2�˵�ȥ·����·
        % ---------------- ɾȥj��ȥ· -------------- %
        for ii=1:2
            j = GT_move{t}(ind,ii); % jΪGT_moveһ���е�ǰ2��
            if j % j=0ʱΪk���֣�j��ȥ·
                % 1 Ǩ�Ƴ�����0
                Fij_c{t}(j,:) = zeros(1,4);
                % 2 divide/split������0
                Fid_c{t}(j,:) = zeros(1,6);
                Fiv_c{t}(j,:) = zeros(1,6);
                % 3 ��ʧ������0
                Fit_c{t}(j) = 0;
                % 4 merge������0
                uu = conflict_pair_last_xy{t}{j};
                for i=1:size(uu,1)
                    Fmj_c{t+1}(uu(i,1),uu(i,2)) = 0;
                end
                disp([ '��', num2str(t), '֡���Ϊ', num2str(j), '����Բȥ·��ɾ����']);
            end
        end
        % ---------------- ɾȥk����· -------------- %
        for ii=3:4
            k = GT_move{t}(ind,ii); % kΪGT_moveһ���еĺ�2��
            if k % k=0ʱΪj��ʧ��k����·
                % 1 Ǩ�������0
                uu = conflict_fij{t}{k};
                for i=1:size(uu,1)
                    Fij_c{t}(uu(i,1),uu(i,2)) = 0;
                end
                % 2 merge�����0
                Fmj_c{t+1}(k,:) = zeros(1,6);
                % 3 divide/split�����0
                uu = conflict_pair_next_xy{t+1}{k};
                for i=1:size(uu,1)
                    Fid_c{t}(uu(i,1),uu(i,2)) = 0;
                    Fiv_c{t}(uu(i,1),uu(i,2)) = 0;
                end
                % 4 ���������0
                Fsj_c{t+1}(k) = 0;
                disp([ '��', num2str(t+1), '֡���Ϊ', num2str(k), '����Բ��·��ɾ����']);
            end
        end
        % ------------------------------------------- %
        
        %% B �����ȷ�Ĺ켣
        % =============== ����µĹ켣 =============== %
        % �ж��¼�
        j1 = rowdata(1);
        j2 = rowdata(2);
        k1 = rowdata(3);
        k2 = rowdata(4);
        
        disp('  �����µĹ켣...');
        
        abcd = num2str(rowdata~=0);
        switch abcd
            case '1  0  0  0'
                Fit_c{t}(j1) = 1;
                disp([ '  ��', num2str(t), '֡��', num2str(j1), '����disappear�¼�']);
            case '0  0  1  0'
                Fsj_c{t+1}(k1) = 1;
                disp([ '  ��', num2str(t+1), '֡��', num2str(k1), '����appear�¼�']);
            case '1  0  1  0'
                mm = find(candidate_fij{t}(j1,:)==k1);
                if isempty(mm)
                    msgbox([ '  ��', num2str(t), '֡��', num2str(j1), 'Ǩ�Ƶ�Ŀ��', num2str(k1), '�����������ڣ�']);
                    candidate_fij{t}(j1,1) = k1; % �޸�����
                end
                Fij_c{t}(j1,mm) = 1;   
                disp([ '  ��', num2str(t), '֡��', num2str(j1), 'Ǩ��Ϊ', '��', num2str(t+1), '֡��', num2str(k1)]);
            case '1  0  1  1'
                % ------------------------------------------------------- %
                % divide/split��Ҫ�ҳ���ϸ��pairλ��,���ȼ���divide������split���ֶ�����
                mmtrue = 0;
                for mm=1:6
                    if isempty( setdiff(candidate_k_next{t}{j1,mm}, [k1 k2]) )
                        mmtrue = mm;
                    end
                end
                % �ж��Ƿ������������ҵ���[k1 k2]
                if mmtrue
%                     if Ellipse{t+1}{k1}.ind_region~=Ellipse{t+1}{k2}.ind_region % ����ͬһ����Ϊ���ѣ�����Ϊ����
                        Fid_c{t}(j1,mmtrue) = 1;
                        disp([ '  ��', num2str(t), '֡��', num2str(j1), '����Ϊ',...
                            '��', num2str(t+1), '֡��', num2str(k1), '��', num2str(k2)]);
%                     else
%                         Fiv_c{t}(j1,mmtrue) = 1;
%                         disp([ '  ��', num2str(t), '֡��', num2str(j1), '����Ϊ',...
%                             '��', num2str(t+1), '֡��', num2str(k1), '��', num2str(k2)]);
%                     end        
                else
                    msgbox(['  GT_move�ĵ�',num2str(t),'֡��',num2str(ind),'�е�divide/split����޷��ں�ѡ��Բ���ҵ���'])
                    candidate_k_next{t}{j1,1} = [k1 k2];
                end
                % ------------------------------------------------------- %
            case '1  1  1  0'
                % ------------------------------------------------------- %
                % merge ��Ҫ���ҳ���Դϸ��pairλ��
                mmtrue = 0;
                for mm=1:6
                    if isempty( setdiff(candidate_k_last{t+1}{k1,mm}, [j1 j2]) )
                        mmtrue = mm;
                    end
                end
                % �ж��Ƿ������������ҵ���[j1 j2]
                if mmtrue
                    Fmj_c{t+1}(k1,mmtrue) = 1;
                    disp([ '  ��', num2str(t), '֡��', num2str(j1), '��', num2str(j2), '�ϲ�Ϊ',...
                        '��', num2str(t+1), '֡��', num2str(k1)]);
                else
                    error(['  GT_move�ĵ�',num2str(t),'֡��',num2str(ind),'�е�merge����޷��ں�ѡ��Բ���ҵ���'])
                end
                % ------------------------------------------------------- %
            otherwise
                msgbox(['  GT_move�ĵ�',num2str(t),'֡��',num2str(ind),'�м�¼���Ǳ�׼��ϸ���¼���']);
                candidate_k_last{t+1}{k1,1} = [j1 j2];
        end
        
        disp(' ');
        % ============================================ %
        
    end
end
     
%% 2��ɾ�����Ϊɾ������Բ����·��ȥ·����Щ�Ǵ����ǣ����޷�����Ϊ��ȷ�ģ�
for t=1:numel(GT_delete)    
    % Ϊ�����������ӿ��ٶ�
    if isempty(GT_delete{t})
        continue;
    end
    
    for ind=1:numel(GT_delete{t})
        j = GT_delete{t}(ind);
        if ~j % j=0������
            continue;
        end

        % ��һ�������ɾ����·
        if t>1
            % 1 Ǩ�������0
            uu = conflict_fij{t-1}{j};
            for i=1:size(uu,1)
                Fij_c{t-1}(uu(i,1),uu(i,2)) = 0;
            end
            % 2 merge�����0
            Fmj_c{t}(j,:) = zeros(1,6);
            % 3 divide/split�����0
            uu = conflict_pair_next_xy{t}{j};
            for i=1:size(uu,1)
                Fid_c{t-1}(uu(i,1),uu(i,2)) = 0;
                Fiv_c{t-1}(uu(i,1),uu(i,2)) = 0;
            end
            % 4 ���������0
            Fsj_c{t}(j) = 0;
        end
        disp([ '��', num2str(t), '֡���Ϊ', num2str(j), '����Բ��·��ɾ����']);

        % �ڶ��������ɾ��ȥ·
        if t<frame
            % 1 Ǩ�Ƴ�����0
            Fij_c{t}(j,:) = zeros(1,4);
            % 2 divide/split������0
            Fid_c{t}(j,:) = zeros(1,6);
            Fiv_c{t}(j,:) = zeros(1,6);
            % 3 ��ʧ������0
            Fit_c{t}(j) = 0;
            % 4 merge������0
            uu = conflict_pair_last_xy{t}{j};
            for i=1:size(uu,1)
                Fmj_c{t+1}(uu(i,1),uu(i,2)) = 0;
            end
        end
        disp([ '��', num2str(t), '֡���Ϊ', num2str(j), '����Բȥ·��ɾ����']);
    end
end
disp(' ');
disp('========================================');
% -----------------------------

%% ����������ɵ�GT
%  Ϊ�˲��ı������㾫�ȵĴ��룬�˴��ѽ������ΪFij����ʽ
Fij = Fij_c;
Fit = Fit_c;
Fid = Fid_c;
Fiv = Fiv_c;
Fmj = Fmj_c;
Fsj = Fsj_c;

% ���������ֶ����divide/split/merge�¼����ٽ��б���
% Ŀǰ�С�����
% 2015.7.22�Ѿ������Զ����ã�����������¼��ֻ���¼split����
% ---- split ---- %
% t  i  j1  j2  mm
% 64 75 


if 1
    disp('�����˹��޸Ĺ����GT���̱���');
    save([ trackpath, '\GT\GT_after_hand_tune\GT_Flow_Variables_New.mat'],...
        'Fij','Fit','Fid','Fiv','Fmj','Fsj');
end

% ������ֲ���4�����6�����ڵ���������ʵ��޸������ڵ�����Ҳ���޸�label2e�����ڴ˱���
if 0
    % ע�⣡���޸���candidate����Ҫ��CX_ILP_Pair_Pre_New��������� conflict_pair��conflict_fij
    tic
    disp('���� predata...');
    save([ trackpath, '\Pair\Pre_data_New.mat'], 'Ellipse', 'candidate_k_last',...
        'candidate_k_next', 'conflict_pair_last_xy', 'conflict_pair_next_xy', 'n', 'num_var', 'num_var_sum',...
        'candidate_fij', 'conflict_fij');
    toc
end

