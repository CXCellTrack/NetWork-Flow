function [ KEXI W OBJ ] = CXSL_Update_W_For_BMRM( A, B, lambda )
%UNTITLED Summary of this function goes here
% ========================================= %
%   CX 2015.6.17
%	���룺 A��B���ͷ���lambda���������� N
%   ���̣� ������ι滮���⣨����14��
%   ����� �ɳڱ���KEXI�����º��Ȩ���� W
% ========================================= %
% 
% ����ά��
dimension = numel(A{1});
% �������
w = sdpvar(dimension, 1, 'full');
kexi = sdpvar(1, 1, 'full'); % ���������kexi�ľ�ֵ���漴��

% Ŀ�꺯����1/N* sum( kexi ) + lambda* ||w||^2/2

obj_fun = kexi + lambda* (w'*w) /2;

% Լ����������Ϊ���ӣ������еĵ�ǰA��B�н���Լ��
t_cur = numel(A);
while isempty( A{t_cur} ) % ���Ŀǰ������Ai����
    t_cur = t_cur -1;
end
    
F = [];
for i=1:t_cur
    F = [ F, A{i}'* w + B{i} <= kexi ];
end
F = [ F, kexi >= 0 ];

% �����ι滮����
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
    




























