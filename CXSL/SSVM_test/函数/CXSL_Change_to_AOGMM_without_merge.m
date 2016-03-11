function CXSL_Change_to_AOGMM_without_merge( dataset, flowvars_path )

% ���������merge��split�ĸ��������merge�����³��֣�split����divide

% ���ҵõ��ĸ��ٽ��ת��Ϊcell track challenge�ϵı�׼��ʽ
% dataset = 'training';
[ segpath, trackpath ] = getpath( dataset );

% ʹ��ȫ�ֱ���
global Fij Fit Fid Fiv Fmj Fsj;
global conflict_fij conflict_pair_last_xy conflict_pair_next_xy;

load( flowvars_path );
load([trackpath, '\Pair\Pre_data_New.mat']);
frame = numel(Fmj);

%% 2����6���¼�����һ�𣬽�merge��split�滻Ϊmove
Final_Table = cell(frame,1); %��Final_Table �Ǽ�¼�������ڵı���Դ���ص����У����Ѷ���Ϊ1������Ϊ0

for t=1:frame
    for j=1:numel(Ellipse{t})
        Final_Table{t}(j,1) = 0; % ����0ռλ�ã������ٽ����޸�
        [ eventIn, eventOut ] = CX_CheckInOut( t, j );
        % ��j�Ľ����ڶ�Ϊ0��˵���Ƿ��õļ�˵����-1���б��
        if isequal(eventIn,zeros(1,6)) && isequal(eventOut,zeros(1,6))
            Final_Table{t}(j,1) = -1;
        end
        % ================== ������ ================== %
        % ================== ������ ================== %
        ev = find(eventOut);
        if isempty(ev)
            continue
        end
        switch ev
            case 1 % Ǩ�Ƴ�ȥ
                nextind = find(Fij{t}(j,:));
                nexte = candidate_fij{t}(j, nextind);
                Final_Table{t}(j,1) = nexte;
            case 2 % ------ ��ʧ ------ %
                Final_Table{t}(j,1) = 0;
            case 3 % ------ ���ѳ�ȥ ------ %
                nextind = find(Fid{t}(j,:));
                sons = candidate_k_next{t}{j, nextind};
                Final_Table{t}(j,1) = sons(1);
                Final_Table{t}(j,2) = sons(2);
                % Ҫ��������ϸ����Դ
                Final_Table{t+1}(sons(1),3) = j;
                Final_Table{t+1}(sons(2),3) = j;
            case 4 % ------ �����ȥ ------ %
                nextind = find(Fiv{t}(j,:));
                sons = candidate_k_next{t}{j, nextind};          
                % �������ǵ�ϸ����˵������merge���������������Ч���޸�Ϊ����
                Final_Table{t}(j,1) = sons(1);
                Final_Table{t}(j,2) = sons(2);
                % Ҫ��������ϸ����Դ
                Final_Table{t+1}(sons(1),3) = j;
                Final_Table{t+1}(sons(2),3) = j; 
            case 5 % ------ merge��ȥ ------ %
                % �ⲿ�ַ�����ڴ������޸�
            case 6 % ------ �³���
                % �Գ�����Ӱ��
        end
    end
end

%% 3������final table����AOG������track_txt
for t=1:frame
    Final_Table{t} = [Final_Table{t}, zeros(size(Final_Table{t},1),3-size(Final_Table{t},2))];
end % �� Final_Table ͳһ��չ��3��
       
disp('���� Final_Table ���� track_txt...');
track_txt = zeros(500,4);
count = 0;
for t=1:frame
    for j=1:size(Final_Table{t},1)  
        tmpt = t;
        tmpj = j;
        if isfield(Ellipse{t}{j}, 'color') || Final_Table{t}(j,1)== -1
            continue
        end
        fprintf('\n���ڸ��� %d:%d -> ', t,j);
        count = count+1;
        Ellipse{t}{j}.color = count; % ����ɫ����ֹ������ظ�
        % -------------------------------------------------------------- %
        lastj = Final_Table{t}(j,3); % ��ϸ��������
        while 1
            if tmpt==frame % �ﵽ���һ֡�����м�¼
                track_txt(count,:) = [count, t-1, tmpt-1, lastj];
                fprintf('=> END ');
                break
            end
            flag = Final_Table{tmpt}(tmpj,1:2)>0;
            if isequal(flag, [1,0])
                tmpj = Final_Table{tmpt}(tmpj,1);
                tmpt = tmpt+1;
                Ellipse{tmpt}{tmpj}.color = count; % ����һ��Ҳ��ɫ
                fprintf('%d:%d -> ', tmpt,tmpj);
            else % ����Ǩ�ƾ�����ʧ�ͷ��ѣ�������ͼ�ߵ���Ҷ�ڵ�
                track_txt(count,:) = [count, t-1, tmpt-1, lastj];
                if isequal(flag, [0,0])
                    fprintf('Death');
                elseif isequal(flag, [1,1])
                    fprintf('Divide => (%d:%d, %d:%d)',...
                            tmpt+1,Final_Table{tmpt}(tmpj,1), tmpt+1,Final_Table{tmpt}(tmpj,2));
                end
                break
            end
        end
        % -------------------------------------------------------------- %
    end
end
fprintf('\n');
ii = 1;
while ii<=size(track_txt,1)
    if isequal(track_txt(ii,:),[0,0,0,0])
        track_txt(ii,:) = [];
    else
        ii = ii +1;
    end
end

for h=1:size(track_txt,1)
    if track_txt(h,4)
        track_txt(h,4) = Ellipse{track_txt(h,2)}{track_txt(h,4)}.color;
    end
end
       
       
% ���� track_txt ���ı�
filetosave = [trackpath(1:end-11), '_RES\res_track.txt'];
fidin = fopen(filetosave,'wt');
if fidin==-1
    warning('RESĿ¼û�У�');
    mkdir([trackpath(1:end-11), '_RES']);
end
for ii=1:size(track_txt,1)
    fprintf(fidin, '%d %d %d %d\n',track_txt(ii,:));
end
fclose(fidin);

%% 3-4��һЩ�������Բ���
% ========================================= %
% man_track000 ����Ϊ3λ���ֲ��У����
gtpath = [trackpath(1:end-11), '_GT\TRA\'];
gt_dir = dir([gtpath,'*.tif']);
for i=1:numel(gt_dir)
    source = [gtpath,gt_dir(i).name];
    desti = [gtpath,'man_track0',gt_dir(i).name(end-5:end)];
    if strcmp(source,desti)
        break
    end
    movefile(source, desti,'f');
end
gtpath = [trackpath(1:end-11), '_GT\SEG\'];
gt_dir = dir([gtpath,'*.tif']);
for i=1:numel(gt_dir)
    source = [gtpath,gt_dir(i).name];
    desti = [gtpath,'man_seg0',gt_dir(i).name(end-5:end)];
    if strcmp(source,desti)
        break
    end
    movefile(source,desti, 'f');
end
% ========================================= %
% % mask000 Ҳһ��
% dirpath = 'E:\datasets\first_edition\training_datasets\N2DL-HeLa\01_RES\';
% gt_dir = dir([dirpath,'*.tif']);
% for i=1:numel(gt_dir)
%     movefile([dirpath,gt_dir(i).name], [dirpath,'mask0',gt_dir(i).name(end-5:end)],'f');
% end 

% % ��gt��Ϊ80֡���Ա�Ƚ�
% man_track = 'E:\datasets\first_edition\training_datasets\N2DL-HeLa\01_GT\TRA\man_track.txt';
% file = load(man_track);
% file(file(:,3)>=79, 3) = 79;
% i = 1;
% while i<=size(file,1)
%     if file(i,2)>79
%         file(i,:) = [];
%     else
%         i = i+1;
%     end
% end
% fidin = fopen(man_track,'wt');
% for ii=1:size(file,1)
%     fprintf(fidin, '%d %d %d %d\n',file(ii,:));
% end

%% 4������Ҷ�ͼ��
saveaddr = [trackpath(1:end-11), '_RES\'];
lunkuopath = [segpath, '\FOI��ȡ����\'];
lunkuodir = dir([lunkuopath,'*.tif']);

for t=1:frame
    disp(['����ͼƬ',num2str(t),'...']);tic
    % ----- ����ԭ�򵥷ָ�ͼ ------- %
    imraw = imread([lunkuopath, lunkuodir(t).name]);
    [height,width] = size(imraw);
    imraw = imfill(imraw,'holes'); % imshow(imraw)
    Label = bwlabeln(imraw);
    % ----------------------------- %
    im = uint16(zeros(height,width));
    
    for j=1:numel(Ellipse{t})  
        if ~isfield(Ellipse{t}{j}, 'color')
            continue
        end
        % �ں�ͼ�ϻ�����Բ
        [im,count] = plot_ellipse_label(im, Label, Ellipse{t}{j}, 1);
        if count==0
            warning([num2str(t),':',num2str(j),' count=0']); 
        end   
    end
    
    name = num2str(t-1);
    if t<=10 % ���ֱ���Ϊ3λ
        name = ['0', name];
    end
    imwrite(im, [saveaddr,'mask0',name,'.tif']);toc
end
        
      







            





