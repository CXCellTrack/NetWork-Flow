function [ KEXI W OBJ ] = CXSL_Update_W_For_BMRM( A, B, lambda )
%UNTITLED Summary of this function goes here
% ========================================= %
%   CX 2015.6.17
%	输入： A、B、惩罚项lambda、样本个数 N
%   过程： 计算二次规划问题（方程14）
%   输出： 松弛变量KEXI、更新后的权向量 W
% ========================================= %
% 
% 特征维数
dimension = numel(A{1});
% 定义变量
w = sdpvar(dimension, 1, 'full');
kexi = sdpvar(1, 1, 'full'); % 这个用所有kexi的均值代替即可

% 目标函数：1/N* sum( kexi ) + lambda* ||w||^2/2

obj_fun = kexi + lambda* (w'*w) /2;

% 约束条件：较为复杂，从所有的当前A和B中进行约束
t_cur = numel(A);
while isempty( A{t_cur} ) % 求出目前包含的Ai个数
    t_cur = t_cur -1;
end
    
F = [];
for i=1:t_cur
    F = [ F, A{i}'* w + B{i} <= kexi ];
end
F = [ F, kexi >= 0 ];

% 求解二次规划问题
options = sdpsettings('verbose',0,'solver','gurobi');
sol = solvesdp( F, obj_fun, options);


if sol.problem == 0
    KEXI = value(kexi);
    W = value(w);
    OBJ = value(obj_fun);
else
    sol.info
    yalmiperror(sol.problem)
end






end
    




























