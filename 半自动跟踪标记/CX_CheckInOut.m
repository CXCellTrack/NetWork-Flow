function CX_CheckInOut( dataset, t, j )

[ ~, trackpath ] = getpath( dataset );

load([ trackpath, '\GT\GT_after_hand_tune\GT_Flow_Variables_New.mat']);
load([ trackpath, '\Pair\Pre_data_New.mat']);