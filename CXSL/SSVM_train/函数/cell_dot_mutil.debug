function sum_f = cell_dot_mutil( feature, z )
%
% ======================================================================= %
%
% 注意，这个函数只用于debug，实际上运行中调用的是写在 CXSL_Calculate_Obj_With_Loss
% 和 CXSL_Calculate_fai_x_zstar 中的内部函数
% 根据已经验证的结果发现：
%       CXSL_Calculate_Obj_With_Loss 中由于数据类型为sdpvar，使用方法二更快
%   即arrayfun运算
%       CXSL_Calculate_fai_x_zstar 中由于数据为double类型，使用for循环反而快
% ======================================================================= %

% feature 为特征向量cell	n(t)*n(t+1)
% z 为流程变量矩阵
% 输出 sum_f 为该cell内所有特征向量之和
% ===================================

% 1. for 循环式计算 4.75秒

% sum_f = zeros( size(feature{1,1}) );
% [h w] = size(feature);
% for i=1:h
%     for j=1:w
%         sum_f = sum_f + feature{i, j}* z(i,j);
%     end
% end

% 2. arrayfun运算 3.10秒

% 这里有2个版本，一个是 CXSL_Calculate_Obj_With_Loss 中计算 fai(x,z^) 的版本
% 由于得到的数据类型为 sdpvar，无法进行cell2mat，因此最后一步只能用sum求和
% 代码如下：
% ================================================= %

sum_f = 0;
ss = numel(feature);
% feature1 = cellfun(@(x)x', feature, 'un',0);
feature2 = reshape(feature, ss, 1);
z1 = reshape(z, ss, 1);
f_z = arrayfun(@(x)feature2{x}*z1(x), 1:ss, 'un',0);
for x=1:ss
    sum_f = sum_f + f_z{x};
end

% ================================================= % 
% 另一个版本是 CXSL_Calculate_fai_x_zstar 中计算 fai(x,z*) 的版本
% 由于数据类型为double，因此可以转换为mat后使用sum求和，代码如下
% 但实际过程中发现这个方法并没有上一个快，因此不采用

% ss = numel(feature);
% feature1 = cellfun(@(x)x', feature, 'un',0);
% feature2 = reshape(feature1, ss, 1);
% z1 = reshape(z, ss, 1);
% f_z = arrayfun(@(x)feature2{x}*z1(x), 1:ss, 'un',0)';
% f_z_mat = cell2mat(f_z);
% sum_f = sum(f_z_mat);

end



