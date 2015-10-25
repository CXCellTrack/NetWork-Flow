function Yhat = assign_y_hat_into_y_i(ev, ind, s_frame, e_frame,Fij,Fit,Fid,Fiv,Fmj,Fsj)

% 根据事件选择Yhat
switch ev 
    case 1
        Yhat = Fij(s_frame(ind):e_frame(ind)-1);
    case 2
        Yhat = Fit(s_frame(ind):e_frame(ind)-1);
    case 3
        Yhat = Fid(s_frame(ind):e_frame(ind)-1);
    case 4
        Yhat = Fiv(s_frame(ind):e_frame(ind)-1);
    case 5
        Yhat = Fmj(s_frame(ind):e_frame(ind)-1);
    case 6
        Yhat = Fsj(s_frame(ind):e_frame(ind)-1);
end
