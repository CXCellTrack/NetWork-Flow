function [ feature feature_p fmin fmax ] = CX_mapminmax( feature, minmaxf )

% ���������� feature ��һ����[-1,1]֮��
% ���ع�һ���� feature����һ�������� feature_p����һ��С��� fmin fmax
% ���Լ���һʱ��Ҫ�õ� fmin fmax
% �Խṹ�� minmax.fij ����ʽ����


if nargin == 1
    % ˵����ѵ��������������Ҫ��һ��
    d_fij = numel(feature{2}{1}); % ����ά��
    tmpmax = zeros(numel(feature), d_fij);
    tmpmin = tmpmax;
    for t=1:numel(feature)
        for d=1:d_fij
            if isempty(feature{t}) % ��һ��Ϊ�������������fmj��fsj
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
        if isempty(feature{t}) % ��һ��Ϊ�������������fmj��fsj
            continue;
        end
        feature{t} = cellfun(@(x)2*(x-fmin)./(fmax-fmin)-1, feature{t}, 'un',0); % ��������һ������-1��1��
        feature_p{t} = cellfun(@(x)[x;1], feature{t}, 'un',0); % ����һ��1�����������
    end
    % ѵ��������һ������
    
elseif nargin == 2
    fmin = 0;
    fmax = 0;
    % ˵���ǲ���������������Ҫ����ѵ����minmax���й�һ��
    feature_p = cell(size(feature));
    for t=1:numel(feature)
        if isempty(feature{t}) % ��һ��Ϊ�������������fmj��fsj
            continue;
        end
        feature{t} = cellfun(@(x)2*(x-minmaxf.min)./(minmaxf.max-minmaxf.min)-1, feature{t}, 'un',0); % ��������һ������-1��1��
        feature_p{t} = cellfun(@(x)[x;1], feature{t}, 'un',0); % ����һ��1�����������
    end
    % ����������һ���
end
    
    
    
    
    
    
    
    
    
    