% 按流程在测试集上从初始分割到生成最终跟踪图片
% 注意，脚本中的参数只能到具体文件中修改

clear;close all;

%% 1、生成椭圆假说 raw_ellipse
run('CX_Network');

%% 2、优化椭圆，生成 pre_data 并保存
run('CX_ILP_Pair_Pre_New');

%% 3、计算特征，生成 feature_plus 并保存
run('CXSL_Combine_feature_all_New');

%% 4、用训练得到的w来求解ILP问题
run('CXSL_ILP_Using_Best_W_New');

%% 5、绘制椭圆的新拟合图
CX_RePlot_Ellipse( 'competition');

%% 6、绘制最终的跟踪结果
run('CX_Visualize_Track_Pair_New');





