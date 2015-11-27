function kernel_ff_all_ev = ssvm_pre_cal_all_kernel_paper(kernel_ff_all_ev, gt_frame, ev, kernel_type_ev, cmd_ev, iind)    

% ͨ����ⷢ�֣���һ��ȷ����֧������phi_i��˵�����Ӧ�ڦ�_i
% ������ṹ����ΪK(phi_i, phi)
% չ��Ϊ��j��k( yj*yk*k(fj,fk) )------��1��

% ����k(fj,fk)����svm�ĺ˺�����ʽ
% ������1��д�ɾ���˷�����ʽΪ��
% 
% Yj'*K(Fj,Fk)*Yk -------------------��2��
% ����YkΪbinvar����
% K(Fj,Fk)ֻ�������йأ�������ΪN*N�֣������л����ظ���
% Yj��2�֣�
%       1���Ǳ�׼��phi_i*��Ӧ��y_i*
%       2����֧������phi_i��Ӧ��y_i
%
% ����Ŀ�꺯����ʱ��2�߶�Ҫ�����
% ��˿�����������е� K(Fj,Fk)*Yk���ٵõ�Yj���ٽ������Ѹ����� K(phi_i, phi)
% �Ӷ�������ѭ���д����ļ���
                
                
%% ����Ԥ��������
[ ~, trackpath ] = getpath( 'training' );
% load([ trackpath, '\Pair\Pre_data_New.mat']); 
load([ trackpath, '\GT\GT_after_hand_tune\GT_Flow_Variables_New.mat']); % �����׼��
load([ trackpath, '\�ṹ��ѧϰ\Feature_Plus_New.mat']); % ��������

sample_feature = cell(1,gt_frame); % ����������Ӧ������

switch ev % �����¼�ѡ��feature�����̱���
    case 1
        disp(['����move�¼��ĺ˺��� ', kernel_type_ev, '...']);
        sample_feature = feature_fij_p(1:gt_frame-1);
    case 2
        disp(['����disappear�¼��ĺ˺��� ', kernel_type_ev, '...']);
        sample_feature = feature_fit_p(1:gt_frame-1);
    case 3
        disp(['����divide�¼��ĺ˺��� ', kernel_type_ev, '...']);
        sample_feature = feature_fid_p(1:gt_frame-1);
    case 4
        disp(['����split�¼��ĺ˺��� ', kernel_type_ev, '...']);
        sample_feature = feature_fiv_p(1:gt_frame-1);
    case 5
        disp(['����merge�¼��ĺ˺��� ', kernel_type_ev, '...']);
        sample_feature = feature_fmj_p(2:gt_frame);
    case 6
        disp(['����appear�¼��ĺ˺��� ', kernel_type_ev, '...']);
        sample_feature = feature_fsj_p(2:gt_frame);
end

%% ���� k(f,f)

% ע�⣬ԭ��y1��y2������ʱ�����Y1'*K(F,F)*Y2 ����<��i,��>
% ����y1��ÿ���и��£���K(F,F)��Y2���ǹ̶���
% ��˿����ֽ� K(F,F) ����������Լ���ѭ����ʱ��

tic; 
fe1 = sample_feature;
fe2 = sample_feature;

s1 = iind(1); % fe1�Ŀ�ʼ֡
e1 = iind(2); % fe1�Ľ���֡
s2 = iind(3); % fe2�Ŀ�ʼ֡
e2 = iind(4); % fe2�Ľ���֡
if any([e1,e2]>gt_frame)
    error('����֡�������֡�ˣ�');
end

for ii=s1:e1-1
    for jj=s2:e2-1 % ii,jj��ʾ��ǰ�ߵ�i֡�ͺ��ߵ�j֡
        if jj<ii
            continue; % ֻ��������������
        end

        disp(['    �����',num2str(ii), '֡���',num2str(jj),'֡�ĺ˺���...']);
        K_mat = sparse(numel(fe1{ii}),numel(fe2{jj})); % ��ϡ������Ÿ�ʡ�ռ�
        
        for hh=1:size(fe1{ii},1) % ������
            if ev==3
                if sum(Fid{ii}(hh,:))==0
                    continue; % ��һ��y*��Ӧ����0��������
                end
                for ss=1:6 % ������
                    kk = (ss-1)*size(fe1{ii},1)+hh; % ������ֱ�����Բ���
                    for mm=1:numel(fe2{jj}) % ����jj�ߵ�mm����Բ
                        % ����svm�ĺ��㷨�������
                        K_mat(kk,mm) = svm_kernel(fe1{ii}{kk}, fe2{jj}{mm}, kernel_type_ev, cmd_ev);% kernel_type_ev, cmd_ev
                    end
                end
            end

        end
        kernel_ff_all_ev{ii,jj} = K_mat;
    end
end
toc
% �����������5��ѭ����������ǰ��֡������֡��ǰ����������������������ǳ���ʱ��





%% �� kernel_result ��չ�����󣨰��Խ��߾���
% ̫���ڴ棬���ǲ������ˣ�ѵ��ʱ��ʱ���Ƽ���
%     for ii=1:gt_frame-1
%         for jj=1:gt_frame-1
%             if isempty(kernel_ff_all{ev}{ii,jj}) && ~isempty(kernel_ff_all{ev}{jj,ii})
%                 kernel_ff_all{ev}{ii,jj} = kernel_ff_all{ev}{jj,ii}';
%             end
%         end
%     end



    
































