function kernel_ff_all_ev = ssvm_pre_cal_all_kernel_paper(kernel_ff_all_ev, gt_frame, ev, kernel_type_ev, cmd_ev, iind)    

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
% load([ trackpath, '\Pair\Pre_data_New.mat']); 
load([ trackpath, '\GT\GT_after_hand_tune\GT_Flow_Variables_New.mat']); % 载入标准答案
load([ trackpath, '\结构化学习\Feature_Plus_New.mat']); % 载入特征

sample_feature = cell(1,gt_frame); % 所有样本对应的特征

switch ev % 根据事件选择feature和流程变量
    case 1
        disp(['计算move事件的核函数 ', kernel_type_ev, '...']);
        sample_feature = feature_fij_p(1:gt_frame-1);
    case 2
        disp(['计算disappear事件的核函数 ', kernel_type_ev, '...']);
        sample_feature = feature_fit_p(1:gt_frame-1);
    case 3
        disp(['计算divide事件的核函数 ', kernel_type_ev, '...']);
        sample_feature = feature_fid_p(1:gt_frame-1);
    case 4
        disp(['计算split事件的核函数 ', kernel_type_ev, '...']);
        sample_feature = feature_fiv_p(1:gt_frame-1);
    case 5
        disp(['计算merge事件的核函数 ', kernel_type_ev, '...']);
        sample_feature = feature_fmj_p(2:gt_frame);
    case 6
        disp(['计算appear事件的核函数 ', kernel_type_ev, '...']);
        sample_feature = feature_fsj_p(2:gt_frame);
end

%% 计算 k(f,f)

% 注意，原本y1，y2都存在时计算的Y1'*K(F,F)*Y2 就是<φi,φ>
% 由于y1在每轮中更新，而K(F,F)，Y2都是固定的
% 因此可以现将 K(F,F) 先算出来，以减少循环的时间

tic; 
fe1 = sample_feature;
fe2 = sample_feature;

s1 = iind(1); % fe1的开始帧
e1 = iind(2); % fe1的结束帧
s2 = iind(3); % fe2的开始帧
e2 = iind(4); % fe2的结束帧
if any([e1,e2]>gt_frame)
    error('结束帧超过最大帧了！');
end

for ii=s1:e1-1
    for jj=s2:e2-1 % ii,jj表示是前者第i帧和后者第j帧
        if jj<ii
            continue; % 只计算右上三角阵！
        end

        disp(['    计算第',num2str(ii), '帧与第',num2str(jj),'帧的核函数...']);
        K_mat = sparse(numel(fe1{ii}),numel(fe2{jj})); % 用稀疏矩阵存放更省空间
        
        for hh=1:size(fe1{ii},1) % 按行来
            if ev==3
                if sum(Fid{ii}(hh,:))==0
                    continue; % 这一行y*对应的是0，则跳过
                end
                for ss=1:6 % 按列算
                    kk = (ss-1)*size(fe1{ii},1)+hh; % 按列拉直后的椭圆序号
                    for mm=1:numel(fe2{jj}) % 后者jj者第mm个椭圆
                        % 调用svm的核算法进行求解
                        K_mat(kk,mm) = svm_kernel(fe1{ii}{kk}, fe2{jj}{mm}, kernel_type_ev, cmd_ev);% kernel_type_ev, cmd_ev
                    end
                end
            end

        end
        kernel_ff_all_ev{ii,jj} = K_mat;
    end
end
toc
% 这里面包含了5重循环：样本、前者帧、后者帧、前者向量、后者向量，计算非常耗时！





%% 将 kernel_result 扩展到满阵（按对角线镜像）
% 太费内存，还是不复制了，训练时临时复制即可
%     for ii=1:gt_frame-1
%         for jj=1:gt_frame-1
%             if isempty(kernel_ff_all{ev}{ii,jj}) && ~isempty(kernel_ff_all{ev}{jj,ii})
%                 kernel_ff_all{ev}{ii,jj} = kernel_ff_all{ev}{jj,ii}';
%             end
%         end
%     end



    
































