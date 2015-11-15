function YKY = ssvm_kernel_train_paper(y_ev, alpha_ev, Kernel_ev, ev, N, ind_train, s_frame, e_frame,...
                                        fij, fit, fid, fiv, fmj, fsj)    

% 设置一个bool变量use_star来指示是y*还是yi                         
use_star = false;                              
if isempty(y_ev)&&isempty(alpha_ev)
    use_star = true;
end
    
[ ~, trackpath ] = getpath( 'training' );
load([ trackpath, '\GT\GT_after_hand_tune\GT_Flow_Variables_New.mat']); % 载入标准答案
y_star = cell(N,1); % y*
switch ev % 根据事件选择feature和流程变量
    case 1
        if use_star
            for ind=1:N
                y_star{ind} = Fij(s_frame(ind):e_frame(ind)-1);
            end
        end
        y_train = fij{ind_train}(s_frame(ind_train):e_frame(ind_train)-1);
    case 2
        if use_star
            for ind=1:N
                y_star{ind} = Fit(s_frame(ind):e_frame(ind)-1);
            end
        end
        y_train = fit{ind_train}(s_frame(ind_train):e_frame(ind_train)-1);
    case 3
        if use_star
            for ind=1:N
                y_star{ind} = Fid(s_frame(ind):e_frame(ind)-1);
            end
        end
        y_train = fid{ind_train}(s_frame(ind_train):e_frame(ind_train)-1);
    case 4
        if use_star
            for ind=1:N
                y_star{ind} = Fiv(s_frame(ind):e_frame(ind)-1);
            end
        end
        y_train = fiv{ind_train}(s_frame(ind_train):e_frame(ind_train)-1);
    case 5
        if use_star
            for ind=1:N
                y_star{ind} = Fmj(s_frame(ind)+1:e_frame(ind));
            end
        end
        y_train = fmj{ind_train}(s_frame(ind_train)+1:e_frame(ind_train)); 
    case 6
        if use_star
            for ind=1:N
                y_star{ind} = Fsj(s_frame(ind)+1:e_frame(ind));
            end
        end
        y_train = fsj{ind_train}(s_frame(ind_train)+1:e_frame(ind_train));
end
    
if use_star
    % 计算ystar*K*y
    YKY = cal_Ystar_K_Y(y_star, y_train, Kernel_ev, ev, N, ind_train, s_frame, e_frame);
else
    % 计算yi*K*y
    YKY = cal_Yhat_K_Y(y_ev, alpha_ev, y_train, Kernel_ev, ev, N, ind_train, s_frame, e_frame);
end




%% 计算 k(phi*,phi) 等同于 Ystar*K*Y
function Y_star_K_Y = cal_Ystar_K_Y(y_star, y_train, Kernel_ev, ev, N, ind_train, s_frame, e_frame)

Y_star_K_Y = 0;
for ind=1:N
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
    
    yKy = binvar(numel(y1), size(thisK,2), 'full');
    for ii=1:numel(y1)
        for jj=1:size(thisK,2)
            tmpy1 = reshape(y1{ii},[],1)'; % 必须都是按列拉直
            tmpy2 = reshape(y2{jj},[],1);
            tmpK = thisK{ii,jj};
            yKy(ii,jj) = tmpy1*tmpK*tmpy2;
        end
    end
    % 这个和还要*该样本中支持向量的个数
    Y_star_K_Y = Y_star_K_Y + sum(yKy(:));
end
% 求<f*, f>无需*alpha，因为alpha和一定为1

%% 计算 k(phi^,phi) 等同于 Yi*K*Y
function Y_i_K_Y = cal_Yhat_K_Y(y_ev, alpha_ev, y_train, Kernel_ev, ev, N, ind_train, s_frame, e_frame)

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

        yKy = binvar(numel(y1), size(thisK,2), 'full');
        for ii=1:numel(y1)
            for jj=1:size(thisK,2)
                tmpy1 = reshape(y1{ii},[],1)'; % 必须都是按列拉直
                tmpy2 = reshape(y2{jj},[],1);
                tmpK = thisK{ii,jj};
                yKy(ii,jj) = tmpy1*tmpK*tmpy2;
            end
        end
        % 这个和还要*该样本中支持向量的对应的aplha
        svK = svK + alpha1*sum(yKy(:));
    end
    
    Y_i_K_Y = Y_i_K_Y + svK;
end
    
































