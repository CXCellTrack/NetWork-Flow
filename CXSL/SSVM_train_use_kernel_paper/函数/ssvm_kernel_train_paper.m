function [ Y_start_K_Y Y_hat_K_Y ] = ssvm_kernel_train_paper(y_ev, alpha_ev, Kernel_ev, ev, N, ind_train, s_frame, e_frame,...
    fij, fit, fid, fiv, fmj, fsj)

    
[ ~, trackpath ] = getpath( 'training' );
load([ trackpath, '\GT\GT_after_hand_tune\GT_Flow_Variables_New.mat']); % 载入标准答案
% global 

y_star = cell(N,1); % y*
switch ev % 根据事件选择feature和流程变量
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
  
% 计算ystar*K*y
Y_start_K_Y = cal_Ystar_K_Y(y_star, y_train, [Kernel_ev], ev, N, ind_train, s_frame, e_frame);
% 计算yi*K*y
Y_hat_K_Y = cal_Yhat_K_Y(y_ev, alpha_ev, y_train, [Kernel_ev], ev, N, ind_train, s_frame, e_frame);


%% 计算 k(phi*,phi) 等同于 Ystar*K*Y
function Y_star_K_Y = cal_Ystar_K_Y(y_star, y_train, Kernel_ev, ev, N, ind_train, s_frame, e_frame)

% global Kernel_ev;

Y_star_K_Y = 0;
for ind=1:N % 对于一个特定的aplha_i
%     disp(['  计算样本',num2str(ind), '...'])
    y1 = y_star{ind};
    y2 = y_train;
    
    % ---------------------------------- %
    s_s = s_frame(ind); % 样本起始帧
    e_s = e_frame(ind)-1; % 样本结束帧
    s_t = s_frame(ind_train); % 训练样本起始帧
    e_t = e_frame(ind_train)-1; % 训练样本结束帧
    % 根据具体事件调整开始和结束帧
    thisK = Kernel_ev(s_s:e_s, s_t:e_t); % 前样本长度*训练样本长度
    for ith=1:size(thisK,1)
        for jth=1:size(thisK,2)
            if isempty(thisK{ith,jth})
                thisK{ith,jth} = Kernel_ev{jth+s_t-1,ith+s_s-1}';% 由于只有右上三角，因此下三角的需要复制得到
            end
        end
    end
    % ---------------------------------- %
    
    % 计算ΣΣystar_i*y*<fi,f>，之后还要乘以对应的αi
    % 原始方法：使用矩阵乘法求和
    yKy = binvar(numel(y1), size(thisK,2), 'full');
    sumyKy = 0;
    for ii=1:numel(y1)
        for jj=1:size(thisK,2) % siz(thisK,2)=numel(y2)
            tmpy1 = reshape(y1{ii},[],1)'; % 必须都是按列拉直
            tmpy2 = reshape(y2{jj},[],1);
            tmpK = thisK{ii,jj};
            % yKy(ii,jj) = tmpy1*tmpK*tmpy2; % 这个太慢
            % 2015.12.31 新方法：无意间发现符号运算时，采用累加比sum(快很多)
            % 猜测可能是因为符号运算的sum很慢的原因
            sumyKy = sumyKy + tmpy1*tmpK*tmpy2;
        end
    end
    % 这个和还要*该样本中支持向量的个数(y* 由于alpha和为1，就不用乘了)
    Y_star_K_Y = Y_star_K_Y + sumyKy; %sum(yKy(:));
    
end
% 求<f*, f>无需*alpha，因为alpha和一定为1

%% 计算 k(phi^,phi) 等同于 Yi*K*Y
function Y_i_K_Y = cal_Yhat_K_Y(y_ev, alpha_ev, y_train, Kernel_ev, ev, N, ind_train, s_frame, e_frame)

% global Kernel_ev;

Y_i_K_Y = 0;
for ind=1:N
    % ---------------------------------- %
    s_s = s_frame(ind); % 样本起始帧
    e_s = e_frame(ind)-1; % 样本结束帧
    s_t = s_frame(ind_train); % 训练样本起始帧
    e_t = e_frame(ind_train)-1; % 训练样本结束帧
    % 根据具体事件调整开始和结束帧
    thisK = Kernel_ev(s_s:e_s, s_t:e_t); % 前样本长度*训练样本长度
    for ith=1:size(thisK,1)
        for jth=1:size(thisK,2)
            if isempty(thisK{ith,jth})
                thisK{ith,jth} = Kernel_ev{jth+s_t-1,ith+s_s-1}';% 由于只有右上三角，因此下三角的需要复制得到
            end
        end
    end
    % ---------------------------------- %
    y2 = y_train;
    
    y_sample1 = y_ev{ind}; % 一个样本中所有支持向量
    alpha_sample1 = alpha_ev{ind}; % 一个样本中所有支持向量对应的alpha
    svK = 0;
    
    for i_in_one=1:numel(y_sample1)
        y1 = y_sample1{i_in_one}; % 该样本中某一个特定的支持向量
        alpha1 = alpha_sample1(i_in_one); % 该样本中某一个特定的支持向量对应的alpha

%         yKy = binvar(numel(y1), size(thisK,2), 'full');
        sumyKy = 0;
        for ii=1:numel(y1)
            for jj=1:size(thisK,2)
                % 采用向量化方法（慢！）
                tmpy1 = reshape(y1{ii},[],1)'; % 必须都是按列拉直
                tmpy2 = reshape(y2{jj},[],1);
                tmpK = thisK{ii,jj};
                % 计算这个最后还需要计算sum(yKy(:))，比for循环慢太多
%                 yKy(ii,jj) = tmpy1*tmpK*tmpy2;
                % 使用累加速度快了20多倍 2015.12.31
                sumyKy = sumyKy + tmpy1*tmpK*tmpy2;
            end
        end
        % 这个和还要*该样本中支持向量的对应的aplha
%         svK = svK + alpha1*sum(yKy(:));
        svK = svK + alpha1*sumyKy;
    end
    
    Y_i_K_Y = Y_i_K_Y + svK;
end
    
































