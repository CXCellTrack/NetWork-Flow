function [ feature feature_p fmin fmax ] = CX_mapminmax( feature, minmaxf )

% 将输入特征 feature 归一化到[-1,1]之间
% 返回归一特征 feature，归一增广特征 feature_p，归一最小最大 fmin fmax
% 测试集归一时需要用到 fmin fmax
% 以结构体 minmax.fij 的形式输入


if nargin == 1
    % 说明是训练特征进来，需要归一化
    d_fij = numel(feature{2}{1}); % 特征维数
    tmpmax = zeros(numel(feature), d_fij);
    tmpmin = tmpmax;
    for t=1:numel(feature)
        for d=1:d_fij
            if isempty(feature{t}) % 第一个为空则跳过：针对fmj和fsj
                continue;
            end
            tmp = cellfun(@(x)x(d), feature{t});
            tmpmax(t,d) = max(tmp(:));
            tmpmin(t,d) = min(tmp(:));
        end
    end
    fmin = (min(tmpmin))';
    fmax = (max(tmpmax))';

    feature_p = cell(size(feature));
    for t=1:numel(feature)
        if isempty(feature{t}) % 第一个为空则跳过：针对fmj和fsj
            continue;
        end
        feature{t} = cellfun(@(x)2*(x-fmin)./(fmax-fmin)-1, feature{t}, 'un',0); % 将特征归一化到【-1，1】
        feature_p{t} = cellfun(@(x)[x;1], feature{t}, 'un',0); % 最后加一个1变成增广特征
    end
    % 训练特征归一化结束
    
elseif nargin == 2
    fmin = 0;
    fmax = 0;
    % 说明是测试特征进来，需要按照训练的minmax进行归一化
    feature_p = cell(size(feature));
    for t=1:numel(feature)
        if isempty(feature{t}) % 第一个为空则跳过：针对fmj和fsj
            continue;
        end
        feature{t} = cellfun(@(x)2*(x-minmaxf.min)./(minmaxf.max-minmaxf.min)-1, feature{t}, 'un',0); % 将特征归一化到【-1，1】
        feature_p{t} = cellfun(@(x)[x;1], feature{t}, 'un',0); % 最后加一个1变成增广特征
    end
    % 测试特征归一完毕
end
    
    
    
    
    
    
    
    
    
    