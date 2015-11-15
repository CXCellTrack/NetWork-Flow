function [ Fkernel_ind ] = calculate_kernel_constraint( islinear, s_frame_ind, e_frame_ind,...
fij_ind, fit_ind, fid_ind, fiv_ind, fmj_ind, fsj_ind )

[ ~, trackpath ] = getpath( 'training' );
load([ trackpath, '\GT\GT_after_hand_tune\GT_Flow_Variables_New.mat']); % 载入标准答案

Fkernel_ind = 0;
ev = 3;
if ~islinear(ev)
    
    switch ev
        case 1
            
        case 2
            
        case 3
            for tt=s_frame_ind:e_frame_ind-1
                for hh=1:size(Fid{tt},1)
                    if sum(Fid{tt}(hh,:))==0 % 如果答案中这一行都为0
                        Fkernel_ind = [ Fkernel_ind, sum(fid_ind{tt}(hh,:))==0 ]; % 则添加约束强制这一行为0
                    end
                end
            end
            
        case 4
            
        case 5
            
        case 6

    end
    

end
    

                        

    
    
    