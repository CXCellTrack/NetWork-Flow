function K = svm_kernel(u, v, kernel_type, cmd)

% �����������svmtrain�ı�﷽����ͬ��
cmdcell = strsplit(' ', cmd);
for ii=1:2:numel(cmdcell)
    if strcmp(cmdcell{ii}, '-d')
        degree = str2num( cmdcell{ii+1} );
    elseif strcmp(cmdcell{ii}, '-g')
        gamma = str2num( cmdcell{ii+1} );
    end
end
        
% �趨����ʽ�Ĵ����͸�˹gamma��ȱʡֵ
if ~exist('degree', 'var')
    degree = 2;
end
if ~exist('gamma', 'var')
    gamma = 1/numel(u);
end

switch kernel_type
    case 'linear'
        K = dot(u,v);
    case 'poly'
        K = (gamma*u'*v)^degree;
    case 'rbf' 
        K = exp(-gamma*norm(u-v)^2);
    case 'sigmoid'
        K = tanh(gamma*u'*v);
    otherwise
        error('��򿪺��� ssvm_kernel �鿴�ɹ�ѡ��ĺ˺�����');
end





