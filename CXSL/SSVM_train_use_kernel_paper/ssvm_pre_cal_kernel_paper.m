function K_Ytrain = ssvm_pre_cal_kernel_paper(N, ev, kernel_type_ev, cmd_ev,...
                    fij, fit, fid, fiv, fmj, fsj, s_frame, e_frame)    

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
load([ trackpath, '\�ṹ��ѧϰ\Feature_New.mat']); % ��������
load([ trackpath, '\GT\GT_after_hand_tune\GT_Flow_Variables_New.mat']); % �����׼��

sample_feature = cell(1,N); % ����������Ӧ������
flowvar = cell(1,N); % ���̱���

switch ev % �����¼�ѡ��feature�����̱���
    case 1
        for ind=1:N
            ss = s_frame(ind);
            ee = e_frame(ind);
            sample_feature{ind} = feature_fij(ss:ee-1);
            flowvar{ind} = fij{ind}(ss:ee-1);
        end
    case 2
        for ind=1:N
            ss = s_frame(ind);
            ee = e_frame(ind);
            sample_feature{ind} = feature_fit(ss:ee-1);
            flowvar{ind} = fit{ind}(ss:ee-1);
        end
    case 3
        for ind=1:N
            ss = s_frame(ind);
            ee = e_frame(ind);
            sample_feature{ind} = feature_fid(ss:ee-1);
            flowvar{ind} = fid{ind}(ss:ee-1);
        end
    case 4
        for ind=1:N
            ss = s_frame(ind);
            ee = e_frame(ind);
            sample_feature{ind} = feature_fiv(ss:ee-1);
            flowvar{ind} = fiv{ind}(ss:ee-1);
        end
    case 5
        for ind=1:N
            ss = s_frame(ind);
            ee = e_frame(ind);
            sample_feature{ind} = feature_fmj(ss+1:ee);
            flowvar{ind} = fmj{ind}(ss+1:ee);
        end
    case 6
        for ind=1:N
            ss = s_frame(ind);
            ee = e_frame(ind);
            sample_feature{ind} = feature_fsj(ss+1:ee);
            flowvar{ind} = fsj{ind}(ss+1:ee);
        end
end
train_feature = sample_feature; % ѵ��ʱѡȡ��������ʵ����������������ѡȡһ��

%% ���� k(f,f)

K_Ytrain = cell(N, N);
tic;
for i_sample=1:N 
    fe1 = sample_feature{i_sample};
%     y1 = y_input{ind};
    for i_train=i_sample:N % ��ΪK(x1,x2)=K(x2,x1),���ֻ������������󼴿�
        disp(['  ��������',num2str(i_sample), '������',num2str(i_train),'�ĺ˺���...']);
        fe2 = train_feature{i_train};
        y2 = flowvar{i_train};
    
        frame_result = cell(numel(fe1), numel(fe2));
        for ii=1:numel(fe1)
            for jj=1:numel(fe2) % ii,jj��ʾ��ǰ�ߵ�i֡�ͺ��ߵ�j֡

                disp(['    ����ǰ�ߵ�',num2str(ii), '֡����ߵ�',num2str(jj),'֡�ĺ˺���...']);
                K_mat = zeros(numel(fe1{ii}),numel(fe2{jj}));
                for kk=1:numel(fe1{ii}) % ǰ��ii֡��kk����Բ
                    for mm=1:numel(fe2{jj}) % ����jj�ߵ�mm����Բ
                        % ����svm�ĺ��㷨�������
                        K_mat(kk,mm) = svm_kernel(fe1{ii}{kk}, fe2{jj}{mm}, kernel_type_ev, cmd_ev);% kernel_type_ev, cmd_ev
                    end
                end
                % ��y1��y2��ֱ
%                 tmpy1 = reshape(y1{ii}, 1, []);
                tmpy2 = reshape(y2{jj}, [], 1);
                % ����൱�� ��i��j( yi*yj*k(fei,fej) )
%                 frame_result(ii,jj) = tmpy1*K_mat*tmpy2;
                frame_result{ii,jj} = {K_mat*tmpy2};
            end
        end
        % ����ͻ�Ҫ*��������֧�������ĸ���
        K_Ytrain{i_sample,i_train} = frame_result;
    end
end
toc
% �����������5��ѭ����������ǰ��֡������֡��ǰ����������������������ǳ���ʱ��

% �� K_Y_train ��չ�����󣨰��Խ��߾���
for iii=1:N
    for jjj=1:N
        if isempty(K_Ytrain{iii,jjj})
            K_Ytrain{iii,jjj} = K_Ytrain{jjj,iii};
        end
    end
end


    
































