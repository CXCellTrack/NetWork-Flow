function press_c_space( hObject, eventdata, handles)
%
%  
global GT_move GT_delete global_x global_y t;

keyPress = get(gcf, 'currentcharacter');

if keyPress=='c'
    % 按下'c'则此bsp记录完毕，光标右移
    msgbox('CSP');
    global_y = global_y + 1;
end

if keyPress=='s'
    % 按下's'则确认当前行为split记录
    msgbox('split');
    GT_move{t}{global_x,5} = 'split';
end

% 第一帧中如果发生了消失事件，则点击椭圆后按下空格键，就会记录
if keyPress==' '
    % 按下空格则此行记录完毕，需要判断事件
    abcd = num2str( ~isemptycell(GT_move{t}(global_x,1:4)) );
    switch abcd
        case '1  0  0  0'
            msgbox('disappear');
        case '0  0  1  0'
            msgbox('appear');
        case '1  0  1  0'
            msgbox('move');
        case '1  0  1  1'
            if isempty(GT_move{t}{global_x,5})
                msgbox('divide');
            else
                msgbox('split');
            end
        case '1  1  1  0'
            msgbox('merge');
        otherwise
            % 错误情况则将此行清空，并重新选取
            GT_move{t}(global_x,:) = cell(1,4);
            global_x = global_x - 1;
            msgbox('error! please re choose'); 
    end
    global_x = global_x+1; % 进入下一行
    global_y = 1; % 光标回复
end

if gcf==1 && strcmp(keyPress, 'd')
    % 如果需要删除该椭圆，则按d，只能在第一帧中处理该事件
    % 标记为删除的椭圆会在后续被移除来路和去路
    GT_delete{t} = [ GT_delete{t}; GT_move{t}(global_x,1)];
    GT_move{t}{global_x,1} = []; % 需要把之前记录的迁移信息销毁
    % 提示进行了delete操作
    msgbox('delete');
end



end

