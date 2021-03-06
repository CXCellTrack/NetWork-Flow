function K = svm_kernel(u, v, kernel_type_ev, cmd_ev)

% 读入参数（与svmtrain的表达方法相同）
cmdcell = strsplit(' ', cmd_ev);
for ii=1:2:numel(cmdcell)
    if strcmp(cmdcell{ii}, '-d')
        degree = str2num( cmdcell{ii+1} );
    elseif strcmp(cmdcell{ii}, '-g')
        gamma = str2num( cmdcell{ii+1} );
    elseif strcmp(cmdcell{ii}, '-s')
        sig = str2num( cmdcell{ii+1} );
    end
end
        
% 设定多项式的次数和高斯gamma的缺省值
if ~exist('degree', 'var')
    degree = 2;
end
if ~exist('gamma', 'var')
    gamma = 1/numel(u);
end

switch kernel_type_ev(1)
    case 'l' % 'linear'
        K = dot(u,v);
    case 'p' % 'poly'
        K = (gamma*u'*v)^degree;
    case 'r' % 'rbf' 
        K = exp(-gamma*norm(u-v)^2);
    case 's' % 'sigmoid'
        K = tanh(sig*u'*v);
    otherwise
        error('请打开函数 ssvm_kernel 查看可供选择的核函数！');
end





