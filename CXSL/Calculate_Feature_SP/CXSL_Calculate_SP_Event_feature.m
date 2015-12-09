% =========================================================================
% 2015.12.9
% 
% ���� Superpixel ��˵�����¼�������
%
% =========================================================================

clear;close all
if 0
    dataset = 'competition';
else
    dataset = 'training';
end
[ segpath trackpath ] = getpath( dataset );

pre_data_addr = [ trackpath, '\Pair\Pre_data_New.mat'];
load( pre_data_addr );

% ԭʼά������
mkdir([trackpath, '\�ṹ��ѧϰ']);
Feature_New_addr = [ trackpath, '\�ṹ��ѧϰ\Feature_New.mat'];
% ��������
Feature_Plus_New_addr = [ trackpath, '\�ṹ��ѧϰ\Feature_Plus_New.mat']; 

%% ���㵥��SP�����������뵽�����struct��
disp('  ���㳬������������...');
tic;
SPF = CXSL_Calculate_SP_feature( SuperPixel, segpath, n ); % SPFΪSuperpixel��feature
toc;

frame = numel(SPF);
frame = 10;

%% ���濪ʼ����ǰ��2֡��Բ��������
disp('  �����ϸ���¼�����...');

%% ================== move ================== Ǩ������ȫ�����������ʽ
%
% һ��12ά����
%
%###########################################  ʱ�仨�� 32�� 
tic;
feature_fij = cell(frame-1,1);
for t=1:frame-1
    for j=1:n(t)
        spj = SPF{t}{j};
        
        for mm=1:4
            % ӳ��Ѱ�� k
            k = candidate_fij{t}(j,mm);
            % ��ѡ��˵spk
            spk = SPF{t+1}{k};
            % =========================================================== %
            diff_size = abs(spk.Area - spj.Area); % �����ľ���ֵ
            distance = norm(spk.Centroid - spj.Centroid); % λ�ò���ʹ��ŷ�Ͼ���
            % һϵ�б�������ʹ�ò����ֵ
            somej = struct2cell(spj);
            somek = struct2cell(spk);
            diff_some = abs(cell2mat(somej(3:11)) - cell2mat(somek(3:11)));
            % ����Ҷ�ֱ��ͼ��EMD
            hist_j = spj.HistIntensity;
            hist_k = spk.HistIntensity;
            % --- �����ⲿ���� emd_hat_gd_metric_mex ���� EM ����
            EMD = emd_hat_gd_metric_mex(hist_j, hist_k, ones(numel(hist_j), numel(hist_j)) ,0 ,1);

            % �ϳ�move�����������һ��1����Ϊ�������� 2015.6.24�������Ƶ�mapmimax��һ���н�������Ĵ���
            feature_fij{t}{j,mm} = [ diff_size; distance; diff_some; EMD ];

        end      
    end
end
toc;

%% ==================  divide ================== 
%
% һ��24ά����
%
%########################################### 
tic;
feature_fid = cell(frame-1,1);
for t=1:frame-1
    for j=1:n(t)
        father = SPF{t}{j};
        
        for mm=1:6
            sons = candidate_k_next{t}{j,mm};
            son1 = SPF{t+1}{sons(1)};
            son2 = SPF{t+1}{sons(2)};
            
            % �ҶȺͲ���
            diff_intensity_sum = father.SumIntensity - son1.SumIntensity - son2.SumIntensity;
            % ����Ͳ���
            diff_size_sum = father.Area - son1.Area - son2.Area;
            % �����
            distance = norm(father.Centroid - son1.Centroid) + norm(father.Centroid - son2.Centroid);
            % �Ƕȣ�ʹ�������ڻ�����Ƕ�
            v1 = son1.Centroid - father.Centroid;
            v2 = son2.Centroid - father.Centroid;
            angle_patten = acosd(dot(v1,v2)/( norm(v1)*norm(v2) ));
            
            % 2����ϸ���Ĳ��죨һϵ�б�����
            some_son1 = struct2cell(son1);
            some_son2 = struct2cell(son2);
            diff_son = abs(cell2mat(some_son1(3:11)) - cell2mat(some_son2(3:11)));
            
            % ���Ĳ�����Ϣ
            some_father = struct2cell(father); % ����2�����������������ӽ���
            father_fe = cell2mat(some_father([1,3:11,13]));
            
            % �ϳɷ��������������һ��1����Ϊ�������� 2015.6.24
            feature_fid{t}{j,mm} = [ diff_intensity_sum; diff_size_sum; distance; angle_patten;...
                diff_son; father_fe ];
        end
    end
end
toc;   

%% ==================  appear ================== 
%
% һ��11ά����
%
%########################################### 
tic
feature_fsj = cell(frame-1,1);
for t=2:frame
    for j=1:n(t)
        spj = SPF{t}{j};
        some = struct2cell(spj); % ����2�����������������ӽ���
        feature = cell2mat(some([1,3:11,13]));
        
        % �����һ��1����Ϊ�������� 2015.6.24
        feature_fsj{t}{j,1} = feature;
    end
end
toc

%% ==================  disappear ================== 
%
% һ��11ά����
%
%########################################### 
tic
feature_fit = cell(frame-1,1);
for t=1:frame-1
    for j=1:n(t)
        spj = SPF{t}{j};
        some = struct2cell(spj); % ����2�����������������ӽ���
        feature = cell2mat(some([1,3:11,13]));
        
         % �����һ��1����Ϊ�������� 2015.6.24
        feature_fit{t}{j,1} = feature;
    end
end
toc

%% ==================  split ================== 
%
% һ��23ά����
%
%########################################### 
tic
feature_fiv = cell(frame-1,1);
for t=1:frame-1
    for j=1:n(t)
        father = SPF{t}{j};
        
        for mm=1:6
            sons = candidate_k_next{t}{j,mm};
            son1 = SPF{t+1}{sons(1)};
            son2 = SPF{t+1}{sons(2)};
            
            % �ҶȺͲ��죨������abs����divide��ͬ��
            diff_intensity_sum = abs(father.SumIntensity - son1.SumIntensity - son2.SumIntensity);
            % ����Ͳ��죨������abs����divide��ͬ��
            diff_size_sum = abs(father.Area - son1.Area - son2.Area);
            % �����
            distance = norm(father.Centroid - son1.Centroid) + norm(father.Centroid - son2.Centroid);
            
            % 2����ϸ���Ĳ��죨һϵ�б�����
            some_son1 = struct2cell(son1);
            some_son2 = struct2cell(son2);
            diff_son = abs(cell2mat(some_son1(3:11)) - cell2mat(some_son2(3:11)));
            
            % ���Ĳ�����Ϣ
            some_father = struct2cell(father); % ����2�����������������ӽ���
            father_fe = cell2mat(some_father([1,3:11,13]));
            
            % �ϳ�split�����������һ��1����Ϊ�������� 2015.6.24
            feature_fiv{t}{j,mm} = [ diff_intensity_sum; diff_size_sum; distance;...
                diff_son; father_fe ];
        end
    end
end
toc

%% ==================  merge ================== 
%
% һ��23ά����
%
%########################################### 
tic
feature_fmj = cell(frame-1,1);
for t=2:frame
    for j=1:n(t)
        father = SPF{t}{j};
        
        for mm=1:6
            sources = candidate_k_last{t}{j,mm};
            son1 = SPF{t-1}{sources(1)};
            son2 = SPF{t-1}{sources(2)};
            
            % �ҶȺͲ��죨������abs����divide��ͬ��
            diff_intensity_sum = abs(father.SumIntensity - son1.SumIntensity - son2.SumIntensity);
            % ����Ͳ��죨������abs����divide��ͬ��
            diff_size_sum = abs(father.Area - son1.Area - son2.Area);
            % �����
            distance = norm(father.Centroid - son1.Centroid) + norm(father.Centroid - son2.Centroid);
            
            % 2����ϸ���Ĳ��죨һϵ�б�����
            some_son1 = struct2cell(son1);
            some_son2 = struct2cell(son2);
            diff_son = abs(cell2mat(some_son1(3:11)) - cell2mat(some_son2(3:11)));
            
            % ���Ĳ�����Ϣ
            some_father = struct2cell(father); % ����2�����������������ӽ���
            father_fe = cell2mat(some_father([1,3:11,13]));
            
            % �ϳ�merge�����������һ��1����Ϊ�������� 2015.6.24
            feature_fmj{t}{j,mm} = [ diff_intensity_sum; diff_size_sum; distance;...
                diff_son; father_fe ];
        end
    end
end
toc

%% ��������һ������-1��1�����䣬�����������(2015.9.29)
tic
disp('����������һ��...');
if strcmp(dataset, 'training')
    [ feature_fij feature_fij_p minmax.fij.min minmax.fij.max ] = CX_mapminmax( feature_fij );
    [ feature_fid feature_fid_p minmax.fid.min minmax.fid.max ] = CX_mapminmax( feature_fid );
    [ feature_fit feature_fit_p minmax.fit.min minmax.fit.max ] = CX_mapminmax( feature_fit );
    [ feature_fiv feature_fiv_p minmax.fiv.min minmax.fiv.max ] = CX_mapminmax( feature_fiv );
    [ feature_fmj feature_fmj_p minmax.fmj.min minmax.fmj.max ] = CX_mapminmax( feature_fmj );
    [ feature_fsj feature_fsj_p minmax.fsj.min minmax.fsj.max ] = CX_mapminmax( feature_fsj );
    save([ trackpath, '\�ṹ��ѧϰ\minmax.mat'], 'minmax');
else
    [ ~, traintrackpath ] = getpath( 'training' ); % ������Ҫ��ѵ���������һ����
    load([ traintrackpath, '\�ṹ��ѧϰ\minmax.mat'], 'minmax');
    [ feature_fij feature_fij_p ] = CX_mapminmax( feature_fij, minmax.fij );
    [ feature_fid feature_fid_p ] = CX_mapminmax( feature_fid, minmax.fid );
    [ feature_fit feature_fit_p ] = CX_mapminmax( feature_fit, minmax.fit );
    [ feature_fiv feature_fiv_p ] = CX_mapminmax( feature_fiv, minmax.fiv );
    [ feature_fmj feature_fmj_p ] = CX_mapminmax( feature_fmj, minmax.fmj );
    [ feature_fsj feature_fsj_p ] = CX_mapminmax( feature_fsj, minmax.fsj );
end
toc

%% ������������
disp('  ��������');
% ԭʼ����ֻ���ڼ��� SVM ��ʱ������
save(Feature_New_addr, 'feature_fid','feature_fij','feature_fit','feature_fiv','feature_fmj','feature_fsj');
% �����������һ��
save(Feature_Plus_New_addr, 'feature_fid_p','feature_fij_p','feature_fit_p','feature_fiv_p','feature_fmj_p','feature_fsj_p');

% ��ά��:    ij   id  sj   it    iv  mj
d_feature = 12 + 24 + 11 + 11 + 23 + 23;

