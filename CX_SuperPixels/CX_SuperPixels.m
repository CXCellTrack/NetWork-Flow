%######################################
%
% 2015.12.16 CX on desk
% 作用：这个脚本用于从初始分割中得到椭圆拟合假说
% 数据存储：椭圆数据保存在 raw_ellipse 中
% 依赖关系： 反复调用 CX_fit 计算每帧图片的拟合结果 
% 2015.6.5（加入了 FOI 提取，使边缘的判断符合标准 ）
%
%######################################
clear;close all;

if 0
    dataset = 'competition'; % 选择训练还是测试
else
    dataset = 'training';
end
[ segpath, ~ ] = getpath( dataset );

raw_path = segpath( 1:max(strfind(segpath,'\'))+2 );
raw_dir = dir([raw_path, '\*.tif']);
bw_dir = dir([ segpath, '\*.tif' ]); % 初始分割tif图片地址

% 输出超像素图片地址 
output_addr1 = [ segpath, '\超像素灰度标签\'];
if ~exist(output_addr1, 'dir')
    mkdir(output_addr1);
end
output_addr2 = [ segpath, '\超像素彩色标签\'];
if ~exist(output_addr2, 'dir')
    mkdir(output_addr2);
end
output_addr3 = [ segpath, '\超像素网格图\'];
if ~exist(output_addr3, 'dir')
    mkdir(output_addr3);
end

% 从外部导入SuperPixel.mat
SP_addr = [ output_addr1, 'SuperPixel.mat' ];
if exist(SP_addr, 'file')
    load(SP_addr);
else
    SuperPixel = cell(length(raw_dir),1);
end
%% 开始循环处理图片 
nsp = 3000; % 200个超像素（PSC为3000个）
rm_small = 20; % 除去小于80像素的label

for frame=151:251 %length(bw_dir)
    tic
    % 生成图片地址
    raw_pic = [raw_path, '\', raw_dir(frame).name];
    bw_pic = [segpath, '\', bw_dir(frame).name];
    disp(['  处理',raw_pic ]);
    % 生成标签图片和超像素cell
    [ SP, new_labels, maskim, RGB_label ] = CX_Generate_SP_in_1( raw_pic, bw_pic, nsp, rm_small );
    toc
    % 保存得到的图片
    disp('  保存结果');
    savename1 = [ output_addr1, bw_dir(frame).name(1:end-4), '_sp.png' ];
    imwrite(uint16(new_labels), savename1);
    savename2 = [ output_addr2, bw_dir(frame).name(1:end-4), '_sp_color.png' ];
    imwrite(RGB_label, savename2);
    savename3 = [ output_addr3, bw_dir(frame).name(1:end-4), '_sp_grid.png' ];
    imwrite(maskim, savename3);
    % 保存超像素假说数据
    SuperPixel{frame} = SP;
    save(SP_addr, 'SuperPixel');

end








