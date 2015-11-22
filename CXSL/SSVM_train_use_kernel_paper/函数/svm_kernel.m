function K = svm_kernel(u, v, kernel_type_ev, cmd_ev)

% �����������svmtrain�ı�﷽����ͬ��
cmdcell = strsplit(' ', cmd_ev);
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

switch kernel_type_ev(1)
    case 'l' % 'linear'
        K = dot(u,v);
    case 'p' % 'poly'
        K = (gamma*u'*v)^degree;
    case 'r' % 'rbf' 
        K = exp(-gamma*norm(u-v)^2);
    case 's' % 'sigmoid'
        K = tanh(gamma*u'*v);
    otherwise
        error('��򿪺��� ssvm_kernel �鿴�ɹ�ѡ��ĺ˺�����');
end





