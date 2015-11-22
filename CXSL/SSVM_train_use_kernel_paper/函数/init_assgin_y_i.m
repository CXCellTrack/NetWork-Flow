function y_i = init_assgin_y_i(trackpath, s_frame, e_frame, ev)

% 为初始的phi分配对应的y（也就是标准答案）

load([ trackpath, '\GT\GT_after_hand_tune\GT_Flow_Variables_New.mat']); % 载入标准答案

switch ev % 根据事件选择feature和流程变量
    case 1
        y_i = Fij(s_frame:e_frame-1);
    case 2
        y_i = Fit(s_frame:e_frame-1);
    case 3
        y_i = Fid(s_frame:e_frame-1);
    case 4
        y_i = Fiv(s_frame:e_frame-1);
    case 5
        y_i = Fmj(s_frame+1:e_frame);
    case 6
        y_i = Fsj(s_frame+1:e_frame);
end