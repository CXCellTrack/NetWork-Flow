function pressBlank( hObject, eventdata, handles)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
global GT_move GT_delete global_x t;

keyPress = get(gcf, 'currentcharacter');


% 第一帧中如果发生了消失事件，则点击椭圆后按下空格键，就会记录
if keyPress==' '

    % 按下空格则此行记录完毕，需要判断事件
    abcd = num2str(GT_move{t}(global_x,:)~=0);
    switch abcd
        case '1  0  0  0'
            msgbox('disappear');
        case '0  0  1  0'
            msgbox('appear');
        case '1  0  1  0'
            msgbox('move');
        case '1  0  1  1'
            msgbox('divide/split');
        case '1  1  1  0'
            msgbox('merge');
        otherwise
            % 错误情况则将此行清空，并重新选取
            GT_move{t}(global_x,:) = [0 0 0 0];
            global_x = global_x - 1;
            msgbox('error! please re choose'); 
    end
    global_x=global_x+1;
end

if gcf==1 && strcmp(keyPress, 'd')
    % 如果需要删除该椭圆，则按d，只能在第一帧中处理该事件
    % 标记为删除的椭圆会在后续被移除来路和去路
    GT_delete{t} = [ GT_delete{t}; GT_move{t}(global_x,1)];
    GT_move{t}(global_x,1) = 0; % 需要把之前记录的迁移信息销毁
    % 提示进行了delete操作
    msgbox('delete');
end



end

