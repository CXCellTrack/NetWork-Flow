function K_Ytrain = ssvm_pre_cal_kernel_paper(N, ev, kernel_type_ev, cmd_ev,...
                    fij, fit, fid, fiv, fmj, fsj, s_frame, e_frame)    

% 通过拆解发现，对一个确定的支持向量phi_i来说，其对应于α_i
% 定义其结构化核为K(phi_i, phi)
% 展开为ΣjΣk( yj*yk*k(fj,fk) )------（1）

% 其中k(fj,fk)就是svm的核函数形式
% 而将（1）写成矩阵乘法的形式为：
% 
% Yj'*K(Fj,Fk)*Yk -------------------（2）
% 其中Yk为binvar变量
% K(Fj,Fk)只与特征有关，因此最多为N*N种（且其中还有重复）
% Yj有2种：
%       1、是标准答案phi_i*对应的y_i*
%       2、是支持向量phi_i对应的y_i
%
% 计算目标函数的时候，2者都要求出来
% 因此可以先求出所有的 K(Fj,Fk)*Yk，再得到Yj后再进行组合迅速求出 K(phi_i, phi)
% 从而避免了循环中大量的计算
                
                
%% 载入预计算数据
[ ~, trackpath ] = getpath( 'training' );
load([ trackpath, '\结构化学习\Feature_New.mat']); % 载入特征
load([ trackpath, '\GT\GT_after_hand_tune\GT_Flow_Variables_New.mat']); % 载入标准答案

sample_feature = cell(1,N); % 所有样本对应的特征
flowvar = cell(1,N); % 流程变量

switch ev % 根据事件选择feature和流程变量
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
train_feature = sample_feature; % 训练时选取的样本，实际是在所有样本中选取一个

%% 计算 k(f,f)

K_Ytrain = cell(N, N);
tic;
for i_sample=1:N 
    fe1 = sample_feature{i_sample};
%     y1 = y_input{ind};
    for i_train=i_sample:N % 因为K(x1,x2)=K(x2,x1),因此只需计算上三角阵即可
        disp(['  计算样本',num2str(i_sample), '与样本',num2str(i_train),'的核函数...']);
        fe2 = train_feature{i_train};
        y2 = flowvar{i_train};
    
        frame_result = cell(numel(fe1), numel(fe2));
        for ii=1:numel(fe1)
            for jj=1:numel(fe2) % ii,jj表示是前者第i帧和后者第j帧

                disp(['    计算前者第',num2str(ii), '帧与后者第',num2str(jj),'帧的核函数...']);
                K_mat = zeros(numel(fe1{ii}),numel(fe2{jj}));
                for kk=1:numel(fe1{ii}) % 前者ii帧第kk个椭圆
                    for mm=1:numel(fe2{jj}) % 后者jj者第mm个椭圆
                        % 调用svm的核算法进行求解
                        K_mat(kk,mm) = svm_kernel(fe1{ii}{kk}, fe2{jj}{mm}, kernel_type_ev, cmd_ev);% kernel_type_ev, cmd_ev
                    end
                end
                % 将y1和y2拉直
%                 tmpy1 = reshape(y1{ii}, 1, []);
                tmpy2 = reshape(y2{jj}, [], 1);
                % 这个相当于 ΣiΣj( yi*yj*k(fei,fej) )
%                 frame_result(ii,jj) = tmpy1*K_mat*tmpy2;
                frame_result{ii,jj} = {K_mat*tmpy2};
            end
        end
        % 这个和还要*该样本中支持向量的个数
        K_Ytrain{i_sample,i_train} = frame_result;
    end
end
toc
% 这里面包含了5重循环：样本、前者帧、后者帧、前者向量、后者向量，计算非常耗时！

% 将 K_Y_train 扩展到满阵（按对角线镜像）
for iii=1:N
    for jjj=1:N
        if isempty(K_Ytrain{iii,jjj})
            K_Ytrain{iii,jjj} = K_Ytrain{jjj,iii};
        end
    end
end


    
































