function [ Y_start_K_Y Y_hat_K_Y ] = ssvm_kernel_train_paper(y_ev, alpha_ev, Kernel_ev, ev, N, ind_train, s_frame, e_frame,...
    fij, fit, fid, fiv, fmj, fsj)

    
[ ~, trackpath ] = getpath( 'training' );
load([ trackpath, '\GT\GT_after_hand_tune\GT_Flow_Variables_New.mat']); % �����׼��
% global 

y_star = cell(N,1); % y*
switch ev % �����¼�ѡ��feature�����̱���
    case 1
        for ind=1:N
            y_star{ind} = Fij(s_frame(ind):e_frame(ind)-1);
        end
        y_train = fij{ind_train}(s_frame(ind_train):e_frame(ind_train)-1);
    case 2
        for ind=1:N
            y_star{ind} = Fit(s_frame(ind):e_frame(ind)-1);
        end
        y_train = fit{ind_train}(s_frame(ind_train):e_frame(ind_train)-1);
    case 3
        for ind=1:N
            y_star{ind} = Fid(s_frame(ind):e_frame(ind)-1);
        end
        y_train = fid{ind_train}(s_frame(ind_train):e_frame(ind_train)-1);
    case 4
        for ind=1:N
            y_star{ind} = Fiv(s_frame(ind):e_frame(ind)-1);
        end
        y_train = fiv{ind_train}(s_frame(ind_train):e_frame(ind_train)-1);
    case 5
        for ind=1:N
            y_star{ind} = Fmj(s_frame(ind)+1:e_frame(ind));
        end
        y_train = fmj{ind_train}(s_frame(ind_train)+1:e_frame(ind_train)); 
    case 6
        for ind=1:N
            y_star{ind} = Fsj(s_frame(ind)+1:e_frame(ind));
        end
        y_train = fsj{ind_train}(s_frame(ind_train)+1:e_frame(ind_train));
end
  
% ����ystar*K*y
Y_start_K_Y = cal_Ystar_K_Y(y_star, y_train, [Kernel_ev], ev, N, ind_train, s_frame, e_frame);
% ����yi*K*y
Y_hat_K_Y = cal_Yhat_K_Y(y_ev, alpha_ev, y_train, [Kernel_ev], ev, N, ind_train, s_frame, e_frame);


%% ���� k(phi*,phi) ��ͬ�� Ystar*K*Y
function Y_star_K_Y = cal_Ystar_K_Y(y_star, y_train, Kernel_ev, ev, N, ind_train, s_frame, e_frame)

% global Kernel_ev;

Y_star_K_Y = 0;
for ind=1:N % ����һ���ض���aplha_i
%     disp(['  ��������',num2str(ind), '...'])
    y1 = y_star{ind};
    y2 = y_train;
    
    % ---------------------------------- %
    s_s = s_frame(ind); % ������ʼ֡
    e_s = e_frame(ind)-1; % ��������֡
    s_t = s_frame(ind_train); % ѵ��������ʼ֡
    e_t = e_frame(ind_train)-1; % ѵ����������֡
    % ���ݾ����¼�������ʼ�ͽ���֡
    thisK = Kernel_ev(s_s:e_s, s_t:e_t); % ǰ��������*ѵ����������
    for ith=1:size(thisK,1)
        for jth=1:size(thisK,2)
            if isempty(thisK{ith,jth})
                thisK{ith,jth} = Kernel_ev{jth+s_t-1,ith+s_s-1}';% ����ֻ���������ǣ���������ǵ���Ҫ���Ƶõ�
            end
        end
    end
    % ---------------------------------- %
    
    % ���㦲��ystar_i*y*<fi,f>��֮��Ҫ���Զ�Ӧ�Ħ�i
    % ԭʼ������ʹ�þ���˷����
    yKy = binvar(numel(y1), size(thisK,2), 'full');
    sumyKy = 0;
    for ii=1:numel(y1)
        for jj=1:size(thisK,2) % siz(thisK,2)=numel(y2)
            tmpy1 = reshape(y1{ii},[],1)'; % ���붼�ǰ�����ֱ
            tmpy2 = reshape(y2{jj},[],1);
            tmpK = thisK{ii,jj};
            % yKy(ii,jj) = tmpy1*tmpK*tmpy2; % ���̫��
            % 2015.12.31 �·���������䷢�ַ�������ʱ�������ۼӱ�sum(��ܶ�)
            % �²��������Ϊ���������sum������ԭ��
            sumyKy = sumyKy + tmpy1*tmpK*tmpy2;
        end
    end
    % ����ͻ�Ҫ*��������֧�������ĸ���(y* ����alpha��Ϊ1���Ͳ��ó���)
    Y_star_K_Y = Y_star_K_Y + sumyKy; %sum(yKy(:));
    
end
% ��<f*, f>����*alpha����Ϊalpha��һ��Ϊ1

%% ���� k(phi^,phi) ��ͬ�� Yi*K*Y
function Y_i_K_Y = cal_Yhat_K_Y(y_ev, alpha_ev, y_train, Kernel_ev, ev, N, ind_train, s_frame, e_frame)

% global Kernel_ev;

Y_i_K_Y = 0;
for ind=1:N
    % ---------------------------------- %
    s_s = s_frame(ind); % ������ʼ֡
    e_s = e_frame(ind)-1; % ��������֡
    s_t = s_frame(ind_train); % ѵ��������ʼ֡
    e_t = e_frame(ind_train)-1; % ѵ����������֡
    % ���ݾ����¼�������ʼ�ͽ���֡
    thisK = Kernel_ev(s_s:e_s, s_t:e_t); % ǰ��������*ѵ����������
    for ith=1:size(thisK,1)
        for jth=1:size(thisK,2)
            if isempty(thisK{ith,jth})
                thisK{ith,jth} = Kernel_ev{jth+s_t-1,ith+s_s-1}';% ����ֻ���������ǣ���������ǵ���Ҫ���Ƶõ�
            end
        end
    end
    % ---------------------------------- %
    y2 = y_train;
    
    y_sample1 = y_ev{ind}; % һ������������֧������
    alpha_sample1 = alpha_ev{ind}; % һ������������֧��������Ӧ��alpha
    svK = 0;
    
    for i_in_one=1:numel(y_sample1)
        y1 = y_sample1{i_in_one}; % ��������ĳһ���ض���֧������
        alpha1 = alpha_sample1(i_in_one); % ��������ĳһ���ض���֧��������Ӧ��alpha

%         yKy = binvar(numel(y1), size(thisK,2), 'full');
        sumyKy = 0;
        for ii=1:numel(y1)
            for jj=1:size(thisK,2)
                % ����������������������
                tmpy1 = reshape(y1{ii},[],1)'; % ���붼�ǰ�����ֱ
                tmpy2 = reshape(y2{jj},[],1);
                tmpK = thisK{ii,jj};
                % ������������Ҫ����sum(yKy(:))����forѭ����̫��
%                 yKy(ii,jj) = tmpy1*tmpK*tmpy2;
                % ʹ���ۼ��ٶȿ���20�౶ 2015.12.31
                sumyKy = sumyKy + tmpy1*tmpK*tmpy2;
            end
        end
        % ����ͻ�Ҫ*��������֧�������Ķ�Ӧ��aplha
%         svK = svK + alpha1*sum(yKy(:));
        svK = svK + alpha1*sumyKy;
    end
    
    Y_i_K_Y = Y_i_K_Y + svK;
end
    
































