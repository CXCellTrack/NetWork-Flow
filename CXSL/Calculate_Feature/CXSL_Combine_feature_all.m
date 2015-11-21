% function CXSL_Combine_feature_all_New( str )
% =========================================================================
% 2015.6.16
% ��������� CXSL_Combine_feature ���������ڣ���Ҫ����ȫ�����̱�����Ӧ������
% CXSL_Combine_feature �Ը��ʽ����˹��ˣ�ֻ�������ٵ��������������ں�����ѧϰ�п��ܳ�����
%
% input: Ellipse  ����Բ����
%        n  �Ǹ�֡�е���Բ��Ŀ     
%
% output: �������¼��������������� \01_2-16_track\�ṹ��ѧϰ\Feature.mat ��
%         'feature_ffd','feature_fid','feature_fij','feature_fit','feature_fiv','feature_fmj','feature_fsj'
%
% =========================================================================
clear;close all
if 1
    dataset = 'competition';
else
    dataset = 'training';
end
[ segpath trackpath ] = getpath( dataset );

pre_data_addr = [ trackpath, '\Pair\Pre_data_New.mat'];

last = max(strfind(trackpath, '\'));
rawpic_addr = trackpath(1:last+2);  % ��Ҫ����ԭʼͼƬ�ĵ�ַ
% ԭʼά������
mkdir([trackpath, '\�ṹ��ѧϰ']);
Feature_New_addr = [ trackpath, '\�ṹ��ѧϰ\Feature_New.mat'];
% ��������
Feature_Plus_New_addr = [ trackpath, '\�ṹ��ѧϰ\Feature_Plus_New.mat']; 

%% ��Ҫ���� CX_ILP_Pair_Pre ����õ���Բ����
load( pre_data_addr );
frame = numel(Ellipse);
% frame = 60;

%% ���㵥����Բ�����������뵽�����struct��
disp('������Բ��������...');
tic;
Ellipse = CXSL_Calculate_Ellipse_feature( Ellipse, rawpic_addr );
toc;

%% ���濪ʼ����ǰ��2֡��Բ��������
disp('�����ϸ���¼�����...');

%% ================== move ================== Ǩ������ȫ�����������ʽ
%
% ����1-4   ��position��size��eccentric��alpha����
% ����5     ��16bin �Ҷ�ֱ��ͼ���� ʹ����EMD���� ������˲�����Ϻ���еĸ��� 2015.6.8��
% dist= emd_hat_gd_metric_mex([0.5;0.3;0.2], [0.1;0.1;0.8] ,ones(3,3) ,0 ,1)
% ����6-8   ���ҶȺͣ��ҶȾ�ֵ���Ҷȱ�׼��Ĳ���
%###########################################  ʱ�仨�� 32�� 
tic;
feature_fij = cell(frame-1,1);
feature_fij_p = cell(frame-1,1);
for t=1:frame-1
    for j=1:n(t)
        
        e_j = Ellipse{t}{j};
        fj_1_4 = struct2array( e_j.feature.geo ); % j������������ʱ����
        for mm=1:4
            % ӳ��Ѱ�� k
            k = candidate_fij{t}(j,mm);
            
            % ��ѡ��Բe_k
            e_k = Ellipse{t+1}{k};
            fk_1_4 = struct2array( e_k.feature.geo ); % k������������ʱ����
            % ����1-4�Ĳ���
            diff_position = norm( fj_1_4(1:2)-fk_1_4(1:2) ); % λ�ò���ʹ��ŷ�Ͼ���
            diff_1_4 = [ diff_position;( fj_1_4(3:end) - fk_1_4(3:end) )']; % size�Ȳ���֮�����
            % ����Ҷ�ֱ��ͼ��EMD
            hist_j = e_j.feature.intensity.hist;
            hist_k = e_k.feature.intensity.hist;
             % --- �����ⲿ���� emd_hat_gd_metric_mex ���� EM ����
            diff_5 = emd_hat_gd_metric_mex(hist_j, hist_k, ones(numel(hist_j), numel(hist_j)) ,0 ,1);
            % ��������6-8�Ĳ���
            fj_6_8 = [ e_j.feature.intensity.sum, e_j.feature.intensity.mean, e_j.feature.intensity.devia ];
            fk_6_8 = [ e_k.feature.intensity.sum, e_k.feature.intensity.mean, e_k.feature.intensity.devia ];
            diff_6_8 = (fj_6_8 - fk_6_8)';
            % �ϳ�move�����������һ��1����Ϊ�������� 2015.6.24�������Ƶ�mapmimax��һ���н�������Ĵ���
            feature_fij{t}{j,mm} = [ diff_1_4; diff_5; diff_6_8 ];

        end      
%         others = setdiff( 1:n(t+1),candidate_k );
%         for ind=1:numel(others)
%             feature_fij{t}{j,others(ind)} = Inf;
%         end   
    end

end
toc;

%% ==================  divide ================== 
%
% ����1-2  ���ҶȾ�ֵ���죬�Ƕ���Ϣ
% ����3-7  ��2��ϸ��ec���졢size���졢�ҶȾ�ֵ�Ĳ��졢shpae_compactness��2ϸ�����֮��/͹������������
% ����8-9   ����ec��ĸϸ��ͨ���������ģ������ҶȾ�ֵ
%########################################### 
tic;
feature_fid = cell(frame-1,1);
feature_fid_p = cell(frame-1,1);
for t=1:frame-1
    for j=1:n(t)
        for mm=1:6
            
            sons = candidate_k_next{t}{j,mm};
            son1 = Ellipse{t+1}{sons(1)};
            son2 = Ellipse{t+1}{sons(2)};
            father = Ellipse{t}{j};
            
            % �ҶȺͲ���
            diff_intensity_sum = father.feature.intensity.mean - son1.feature.intensity.sum - son2.feature.intensity.sum;
            % �Ƕȣ�ʹ�������ڻ�����Ƕ�
            v1 = son1.feature.geo.position - father.feature.geo.position;
            v2 = son2.feature.geo.position - father.feature.geo.position;
            angle_patten = acosd(dot(v1,v2)/( norm(v1)*norm(v2) ));
            
            % 2����ϸ����ec����
            diff_sons_ec = abs( son1.feature.geo.eccentric - son2.feature.geo.eccentric );
            % 2����ϸ����size������Ͻ�С���Ǹ�size
            diff_sons_size = abs( son1.feature.geo.size - son2.feature.geo.size )/...
                min(son1.feature.geo.size, son2.feature.geo.size);
            % 2����ϸ���ĻҶȾ�ֵ����
            diff_sons_intensity = abs( son1.feature.intensity.mean - son2.feature.intensity.mean );
            % ϸ����shape compactnessָ��
            shape_compact = pi*(son1.a*son1.b + son2.a*son2.b) /CX_Convhull(son1, son2);
            
            % ��ec
            father_ec = father.feature.geo.eccentric;
            % ���ҶȾ�ֵ
            father_intensity = father.feature.intensity.mean;
            % �ϳɷ��������������һ��1����Ϊ�������� 2015.6.24
            feature_fid{t}{j,mm} = [ diff_intensity_sum, angle_patten, diff_sons_ec, diff_sons_size,...
                diff_sons_intensity, shape_compact, father_ec, father_intensity ]';
        end
    end
end
toc;   

%% ==================  appear ================== 
%
% ����1-2  ��size��dist2border
% ����3-5  ���ҶȺͣ��ҶȾ�ֵ���Ҷȱ�׼��
%
%########################################### 
tic
feature_fsj = cell(frame-1,1);
feature_fsj_p = cell(frame-1,1);
for t=2:frame
    for j=1:n(t)
        e_j = Ellipse{t}{j};
        s_d = [ e_j.feature.geo.size, e_j.feature.dist2border ];
        feature_intensity = [ e_j.feature.intensity.sum, e_j.feature.intensity.mean, e_j.feature.intensity.devia ];
        % �����һ��1����Ϊ�������� 2015.6.24
        feature_fsj{t}{j,1} = [ s_d, feature_intensity ]';
    end
end
toc

%% ==================  disappear ================== 
%
% ����1-2  ��size��dist2border
% ����3-5  ���ҶȺͣ��ҶȾ�ֵ���Ҷȱ�׼��
%
%########################################### 
tic
feature_fit = cell(frame-1,1);
feature_fit_p = cell(frame-1,1);
for t=1:frame-1
    for j=1:n(t)
        e_j = Ellipse{t}{j};
        s_d = [ e_j.feature.geo.size, e_j.feature.dist2border ];
        feature_intensity = [ e_j.feature.intensity.sum, e_j.feature.intensity.mean, e_j.feature.intensity.devia ];
         % �����һ��1����Ϊ�������� 2015.6.24
        feature_fit{t}{j,1} = [ s_d, feature_intensity ]';
    end
end
toc

%% ==================  split ================== 
%
% ����1-2  ���ҶȺͲ��졢size�Ͳ���
% ����3-4  ����ϸ��shape_compactness����ϸ��size����
% ����5    ��Դϸ������϶�ָ��hd��ɾȥ��2015.9.29
%########################################### 
tic
feature_fiv = cell(frame-1,1);
feature_fiv_p = cell(frame-1,1);
for t=1:frame-1
    for j=1:n(t)
        for mm=1:6
            
            couples = candidate_k_next{t}{j,mm};
            brother1 = Ellipse{t+1}{couples(1)};
            brother2 = Ellipse{t+1}{couples(2)};
            father = Ellipse{t}{j};
            
            % �ҶȺͲ���
            diff_intensity_sum = father.feature.intensity.sum - brother1.feature.intensity.sum - brother2.feature.intensity.sum;
            % size�Ͳ���
            diff_size_sum = father.feature.geo.size - brother1.feature.geo.size - brother2.feature.geo.size;
            
            % ϸ����shape compactnessָ��
            shape_compact = pi*(brother1.a*brother1.b + brother2.a*brother2.b) / CX_Convhull(brother1, brother2);
            % ��ϸ��size����
            diff_sons_size = abs( brother1.feature.geo.size - brother2.feature.geo.size );
            
            % Դϸ������϶�ָ��hd
%             source_fitness = father.hd;
            % �����һ��1����Ϊ�������� 2015.6.24
            feature_fiv{t}{j,mm} = [ diff_intensity_sum, diff_size_sum, shape_compact, diff_sons_size ]';
        end
    end
end
toc

%% ==================  merge ================== 
%
% ����1-2  ���ҶȺͲ��졢size�Ͳ���
% ����3-4  ����ϸ��shape_compactness��Դϸ��size����
% ����5    ����ϸ������϶�ָ��hd��ɾȥ��2015.9.29
%########################################### 
tic
feature_fmj = cell(frame-1,1);
feature_fmj_p = cell(frame-1,1);
for t=2:frame
    for j=1:n(t)
        for mm=1:6

            sources = candidate_k_last{t}{j,mm};
            source1 = Ellipse{t-1}{sources(1)};
            source2 = Ellipse{t-1}{sources(2)};
            father = Ellipse{t}{j};
            
            % �ҶȺͲ���
            diff_intensity_sum = father.feature.intensity.sum- source1.feature.intensity.sum - source2.feature.intensity.sum;
            % size�Ͳ���
            diff_size_sum = father.feature.geo.size - source1.feature.geo.size - source2.feature.geo.size;
            
            % ϸ����shape compactnessָ��
            shape_compact = pi*(source1.a*source1.b + source2.a*source2.b) / CX_Convhull(source1, source2);
            % Դϸ��size����
            diff_sons_size = abs( source1.feature.geo.size - source2.feature.geo.size );
            
            % ��ϸ������϶�ָ��hd
%             source_fitness = father.hd;
            % �����һ��1����Ϊ�������� 2015.6.24
            feature_fmj{t}{j,mm} = [ diff_intensity_sum, diff_size_sum, shape_compact, diff_sons_size ]';
        end
    end
end
toc

%% ==================  false detetion ================== 
%
% ����1     ��size
% ����2-4   ���ҶȺͣ��ҶȾ�ֵ���Ҷȱ�׼��, dist2border
%###########################################
feature_ffd = cell(frame-1,1);
feature_ffd_p = cell(frame-1,1);
for t=1:frame-1
    for j=1:n(t)
        p_s = Ellipse{t}{j}.feature.geo.size;
        feature_intensity = [ Ellipse{t}{j}.feature.intensity.sum, Ellipse{t}{j}.feature.intensity.mean, Ellipse{t}{j}.feature.intensity.devia ];
        % �����һ��1����Ϊ�������� 2015.6.24
        feature_ffd{t}{j,1} = [ p_s, feature_intensity, Ellipse{t}{j}.feature.dist2border ]';
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

     
        
 
%% ����������ڼ���2����Բ����С͹������ε����
% ��Ҫ�������ú��� convhull ���м���
if 0 
% function Area = CX_Convhull(e1, e2)
     
% �ȼ���e1��Բ�ܵ㣬ȡ361����
alpha1 = e1.alpha;
a = e1.a;
b = e1.b;
x0 = e1.x0;
y0 = e1.y0;

c=cosd(alpha1);
s=sind(alpha1);
polar_angle=linspace(0,360,50);
xq= a*cosd(polar_angle);
yq= b*sind(polar_angle);
xn1=xq*c-yq*s+x0;
yn1=xq*s+yq*c+y0;

% �ټ���e2��Բ�ܵ㣬ȡ20����
alpha1 = e2.alpha;
a = e2.a;
b = e2.b;
x0 = e2.x0;
y0 = e2.y0;

c=cosd(alpha1);
s=sind(alpha1);
polar_angle=linspace(0,360,50);
xq= a*cosd(polar_angle);
yq= b*sind(polar_angle);
xn2=xq*c-yq*s+x0;
yn2=xq*s+yq*c+y0;

xx = [xn1 xn2];
yy = [yn1 yn2];

[~, Area] = convhull(xx, yy);       
        
end

