function CXSL_Change_to_AOGMM( flowvars_path )


% ���ҵõ��ĸ��ٽ��ת��Ϊcell track challenge�ϵı�׼��ʽ
dataset = 'training';
[ ~, trackpath ] = getpath( dataset );
% result_addr = [trackpath, '\GT\GT_after_hand_tune\'];

% ʹ��ȫ�ֱ���
global Fij Fit Fid Fiv Fmj Fsj;
global conflict_fij conflict_pair_last_xy conflict_pair_next_xy n;
% fig_dir = dir([fig_addr, '*.fig']); // ����Ҫ�õ����Ƴ���fig
% load([result_addr, 'GT_Flow_Variables_New.mat']);
% load([trackpath, '\�²��Խ����¼\local_b.mat']);

load( flowvars_path );
load([trackpath, '\Pair\Pre_data_New.mat']);
frame = numel(Fmj);

%% 1���Ƚ�split��mergeת��Ϊmove

for t=1:frame
    [x,ii] = find(Fmj{t});
    if isempty(x)
        continue
    end
    for i_n=1:numel(x) % ��x��ind�е�ÿһ�Խ��в���
        j = x(i_n); tmpii = ii(i_n);
        fprintf(['\n', num2str(t),'��',num2str(j), '��merge���ģ�']);
        source = candidate_k_last{t}{j, tmpii};
        % ���ô���Բ���ӱ�ǣ�˵������Բ��2����ɵ�
        ss1 = Ellipse{t-1}{source(1)};
        ss2 = Ellipse{t-1}{source(2)};
        nowcount = numel(Ellipse{t});
        Ellipse{t}{nowcount+1} = ss1; Ellipse{t}{nowcount+1}.ext = 1; % ����ext��ǣ�˵�����¼ӵ�
        Ellipse{t}{nowcount+2} = ss2; Ellipse{t}{nowcount+2}.ext = 1;
        Ellipse{t}{j}.twocells = [nowcount+1, nowcount+2];        
        % ������
        tmpt = t;
        tmpj = j;
        while 1
            [ ~, eventOut ] = CX_CheckInOut( tmpt, tmpj );
            ev = find(eventOut);
            if isempty(ev)
                break
            end
            if ev==1 % �����Ǩ�ƵĻ�����Ҫ����һ������ԲҲ���ϱ��
                nextind = find(Fij{tmpt}(tmpj,:));
                nexte = candidate_fij{tmpt}(tmpj, nextind);
                fprintf(' -> %d: %d',tmpt+1,nexte);
                nowcount = numel(Ellipse{tmpt+1});
                Ellipse{tmpt+1}{nowcount+1} = ss1; Ellipse{tmpt+1}{nowcount+1}.ext = 1; % ����ext��ǣ�˵�����¼ӵ�
                Ellipse{tmpt+1}{nowcount+2} = ss2; Ellipse{tmpt+1}{nowcount+2}.ext = 1; % ����ext��ǣ�˵�����¼ӵ�
                Ellipse{tmpt+1}{nexte}.twocells = [nowcount+1, nowcount+2];            
                tmpt = tmpt + 1; % ����t��tmpx
                tmpj = nexte;
            else
                break
            end
        end          
    end
end

%% 2����6���¼�����һ�𣬽�merge��split�滻Ϊmove
Final_Table = cell(frame,1); %��Final_Table �Ǽ�¼�������ڵı���Դ���ص����У����Ѷ���Ϊ1������Ϊ0

for t=1:frame
    for j=1:numel(Ellipse{t})
        if isfield(Ellipse{t}{j}, 'ext') % ��ext��ǣ�˵�����¼ӵģ��޷�check in out
            continue
        end
        [ eventIn, eventOut ] = CX_CheckInOut( t, j );
        % ��j�Ľ����ڶ�Ϊ0��˵���Ƿ��õļ�˵����-1���б��
        if isequal(eventIn,zeros(1,6)) && isequal(eventOut,zeros(1,6))
            Final_Table{t}(j,1) = -1;
        end
        % ================== ������ ================== %
        ev = find(eventIn);
        if ~isempty(ev)
            switch ev
                case 5 % merge���������޸�Ϊmove
                    ind = find(Fmj{t}(j,:));
                    % ע��FT��Ҫ��-1���������ʹ��
                    Final_Table{t}(j,1) = -1;
                    % ---------------------------------- %
                    source = candidate_k_last{t}{j, ind};
                    s1 = Ellipse{t-1}{source(1)};
                    s2 = Ellipse{t-1}{source(2)};
                    if ~isfield(Ellipse{t}{j}, 'twocells')
                        error([num2str(t), '-', num2str(j), '��merge����ȴû��twocells��ǣ�']);
                    end
                    e1 = Ellipse{t}{Ellipse{t}{j}.twocells(1)};
                    e2 = Ellipse{t}{Ellipse{t}{j}.twocells(2)};
                    method = cal_dist(s1,s2,e1,e2);
                    if method
                        Final_Table{t-1}(source(1),1) = Ellipse{t}{j}.twocells(1);
                        Final_Table{t-1}(source(2),1) = Ellipse{t}{j}.twocells(2);
                    else
                        Final_Table{t-1}(source(2),1) = Ellipse{t}{j}.twocells(1);
                        Final_Table{t-1}(source(1),1) = Ellipse{t}{j}.twocells(2);
                    end
            end
        end
        % ================== ������ ================== %
        ev = find(eventOut);
        if isempty(ev)
            continue
        end
        switch ev
            case 1 % Ǩ�Ƴ�ȥ����Ҫ��2�����
                nextind = find(Fij{t}(j,:));
                nexte = candidate_fij{t}(j, nextind);
                if isfield(Ellipse{t}{j}, 'twocells') && isfield(Ellipse{t+1}{nexte}, 'twocells') % ��2����Ϊ��ϸ��
                    % ע��FT��Ҫ��-1���������ʹ��
                    Final_Table{t}(j,1) = -1;
                    Final_Table{t+1}(nexte,1) = -1;
                    % ---------------------------------- %
                    s1 = Ellipse{t}{Ellipse{t}{j}.twocells(1)};
                    s2 = Ellipse{t}{Ellipse{t}{j}.twocells(2)};
                    e1 = Ellipse{t+1}{Ellipse{t+1}{nexte}.twocells(1)};
                    e2 = Ellipse{t+1}{Ellipse{t+1}{nexte}.twocells(2)};
                    method = cal_dist(s1,s2,e1,e2);
                    if method
                        Final_Table{t}(Ellipse{t}{j}.twocells(1),1) = Ellipse{t+1}{nexte}.twocells(1);
                        Final_Table{t}(Ellipse{t}{j}.twocells(2),1) = Ellipse{t+1}{nexte}.twocells(2);
                    else
                        Final_Table{t}(Ellipse{t}{j}.twocells(2),1) = Ellipse{t+1}{nexte}.twocells(1);
                        Final_Table{t}(Ellipse{t}{j}.twocells(1),1) = Ellipse{t+1}{nexte}.twocells(2);
                    end
                else
                    % ��û�д�ϸ��
                    Final_Table{t}(j,1) = nexte;
                end
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
                if isfield(Ellipse{t}{j}, 'twocells') % ��������Ǵ�ϸ�������޸�ΪǨ��
                    % ע��FT��Ҫ��-1���������ʹ��
                    Final_Table{t}(j,1) = -1;
                    % ---------------------------------- %
                    s1 = Ellipse{t}{Ellipse{t}{j}.twocells(1)};
                    s2 = Ellipse{t}{Ellipse{t}{j}.twocells(2)};
                    e1 = Ellipse{t+1}{sons(1)};
                    e2 = Ellipse{t+1}{sons(2)};
                    method = cal_dist(s1,s2,e1,e2);
                    if method
                        Final_Table{t}(Ellipse{t}{j}.twocells(1),1) = sons(1);
                        Final_Table{t}(Ellipse{t}{j}.twocells(2),1) = sons(2);
                    else
                        Final_Table{t}(Ellipse{t}{j}.twocells(2),1) = sons(1);
                        Final_Table{t}(Ellipse{t}{j}.twocells(1),1) = sons(2);
                    end
                else % �������ǵ�ϸ����˵������merge���������������Ч���޸�Ϊ����
                    Final_Table{t}(j,1) = sons(1);
                    Final_Table{t}(j,2) = sons(2);
                    % Ҫ��������ϸ����Դ
                    Final_Table{t+1}(sons(1),3) = j;
                    Final_Table{t+1}(sons(2),3) = j;
                end
            case 5 % ------ merge��ȥ ------ %
                % �ⲿ�ַ�����ڴ������޸�
            case 6 % ------ �³���
                % �Գ�����Ӱ��
        end
    end
end

%% ֮ǰ��˼·�������ã�

% % ����ȥ��˼·���Ȼ�ͼ����tѭ��������ƻҶ���Բ�����θ�����ɫ
% % ����Ellipse�д��ϻҶȱ�ǩ���е�������visualize����
% % ������ͼ�����������޻�ͼ
% count = 0;
% track_txt = zeros(500,4);
% for t=1:frame
%     for j=1:numel(Ellipse{t})
%         if ~isfield(Ellipse{t}{j}, 'color') % ֻ���³��ֵ�ϸ��������
%             count = count+1;
%             Ellipse{t}{j}.color = count;
%             tmpt = t;
%             tmpj = j;
%             while tmpt<frame
%                 [ eventIn, eventOut ] = CX_CheckInOut( tmpt, tmpj );
%                 % ------------------------------------------------------- %
%                 if tmpt==t % ���ԭϸ�������
%                     ev = find(eventIn);
%                     father = 0;
%                     switch ev
%                         case 1: % Ǩ������
%                             error([num2str(t), '-', num2str(j), '����Ǩ������˵������������©��']);
%                         case 3: % ���Ѷ��������ҵ�ĸϸ��
%                             [f_j, f_ii] = find(Fid{t-1});
%                             for ii=1:numel(f_j)
%                                 if any(candidate_k_next(f_j(ii), f_ii(ii))==j)
%                                     father = f_j(ii);
%                                 end
%                             end
%                             if ~father
%                                 error([num2str(t), '-', num2str(j), '�Ƿ��Ѷ�����ȴ�Ҳ���ĸϸ����']);
%                             end
%                             if ~isfield(Ellipse{t-1}{father}, 'color')
%                                 error([num2str(t-1), '-', num2str(father), '��ĸϸ��ȴû�б���ɫ��']);
%                             end
%                             father_color = Ellipse{t-1}{father}.color; % �����������ϸ�������һλ
%                     end
%                 end
%                 % ------------------------------------------------------- %
%                 % ������
%                 ev = find(eventOut);
%                 if t~=frame
%                     assert(numel(ev)==1)
%                 end
%                 switch ev
%                     case 1: % �����Ǩ�ƵĻ�����Ҫ����һ������ԲҲ���ϱ��
%                         nextind = find(Fij{tmpt}(tmpj,:));
%                         nexte = candidate_fij{tmpt}(tmpj,nextind);
%                         Ellipse{t+1}{nexte}.color = count;
%                         tmpt = tmpt + 1; % ����t��tmpx
%                         tmpj = nexte;
%                     case 2: % ������ʧ
%                         track_txt(count,:) = [count, t-1, tmpt-1, father_color];
%                     case 3: % ���з���
%                         track_txt(count,:) = [count, t-1, tmpt-1, father_color];
%                     case 4: % ���з���
%                     case 5: % ����merge

%% 3������final table����AOG������track_txt
for t=1:frame
    Final_Table{t} = [Final_Table{t}, zeros(size(Final_Table{t},1),3-size(Final_Table{t},2))];
end % �� Final_Table ͳһ��չ��3��
       
% file = 'C:\Users\Administrator\Desktop\tracklets.txt'; % ÿһ���켣�ĸ���·��
% delete(file);
% diary(file);
% diary on

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
%         fprintf('\n���ڸ��� %d:%d -> ', t,j);
        count = count+1;
        Ellipse{t}{j}.color = count; % ����ɫ����ֹ������ظ�
        % -------------------------------------------------------------- %
        lastj = Final_Table{t}(j,3); % ��ϸ��������
        while 1
            if tmpt==frame % �ﵽ���һ֡�����м�¼
                track_txt(count,:) = [count, t-1, tmpt-1, lastj];
%                 fprintf('=> END ');
                break
            end
            flag = Final_Table{tmpt}(tmpj,1:2)>0;
            if isequal(flag, [1,0])
                tmpj = Final_Table{tmpt}(tmpj,1);
                tmpt = tmpt+1;
                Ellipse{tmpt}{tmpj}.color = count; % ����һ��Ҳ��ɫ
%                 fprintf('%d:%d -> ', tmpt,tmpj);
            else % ����Ǩ�ƾ�����ʧ�ͷ��ѣ�������ͼ�ߵ���Ҷ�ڵ�
                track_txt(count,:) = [count, t-1, tmpt-1, lastj];
                if isequal(flag, [0,0])
%                     fprintf('Death');
                elseif isequal(flag, [1,1])
%                     fprintf('Divide => (%d:%d, %d:%d)',...
%                         tmpt+1,Final_Table{tmpt}(tmpj,1), tmpt+1,Final_Table{tmpt}(tmpj,2));
                end
                break
            end
        end
        % -------------------------------------------------------------- %
    end
end

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
fprintf('\n');
% diary off             
       
% ���� track_txt ���ı�
filetosave = [trackpath(1:end-11), '_RES\res_track.txt'];
fidin = fopen(filetosave,'wt');
for ii=1:size(track_txt,1)
    fprintf(fidin, '%d %d %d %d\n',track_txt(ii,:));
end
fclose(fidin);

%% 4������Ҷ�ͼ��
% % man_track000 ����Ϊ3λ���ֲ��У����
% dirpath = 'E:\datasets\first_edition\training_datasets\N2DL-HeLa\01_GT\TRA\';
% gt_dir = dir([dirpath,'*.tif']);
% for i=1:numel(gt_dir)
%     movefile([dirpath,gt_dir(i).name], [dirpath,'man_track0',gt_dir(i).name(end-5:end)],'f');
% end
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

height = 700;
width = 1100;
savedir = [trackpath(1:end-11), '_RES\'];

for t=1:frame
    disp(['����ͼƬ',num2str(t),'...']);tic
    im = uint16(zeros(height,width));
    for j=1:numel(Ellipse{t})  
        if ~isfield(Ellipse{t}{j}, 'color')
            continue
        end
        [im,count] = plot_ellipse_label(im, Ellipse{t}{j});
        if count==0
            warning([num2str(t),':',num2str(j),' count=0']); 
        end   
    end
    
    name = num2str(t-1);
    if t<=10 % ���ֱ���Ϊ3λ
        name = ['0', name];
    end
    imwrite(im, [savedir,'mask0',name,'.tif']);toc
end
        
      







            





