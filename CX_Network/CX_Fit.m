%######################################
%
% 2015.6.5 CX on desk
% ���ã�������������ڱ� CX_Network ����
% 2015.6.5�������� FOI ��ȡ��ʹ��Ե���жϷ��ϱ�׼ ��
%
%######################################

function [ e ,lunkuo ]= CX_Fit( frame, pic ,iteration_num, tolerance, remove_small )   %%iteration_num�������˲�����
    
tic;

    %% 1.����ͼƬ
%     error(nargchk(4, 4, nargin));
    if ischar(pic)
        fprintf('���ڴ��� %s ...�˲����� %d\n',pic, iteration_num);  
        im = imread(pic);
    else
        im = pic;
    end
    
    %% 2.��ֵ��
    bw = im2bw(im);
    
    %% 3.����foi(����Ȥ������),�����ֻ�ȫ���������ڵ����£������ɾȥ 2015.6.5 ����welcomeҪ���ƶ�
%   ? E = 50 for the datasets DIC-C2DH-HeLa, Fluo-C2DL-MSC, Fluo-C3DH-H157, Fluo-N2DHGOWT1,
%   Fluo-N3DH-CE, Fluo-N3DH-CHO, and PhC-C2DH-U373,
%   ? E = 25 for the datasets Fluo-C3DL-MDA231, Fluo-N2DL-HeLa, and PhC-C2DL-PSC,
%   ? E = 0 for the simulated datasets (i.e., Fluo-N2DH-SIM, Fluo-N2DH-SIM+, Fluo-N3DH-SIM, and
%   Fluo-N3DH-SIM+) and for the dataset Fluo-N3DL-DRO.
    foi = 0;
    % �����ڲ�����
    bw = Remove_Cells_Out_Of_FOI( bw, foi );

    %% 4.������̬ѧ�ϵ�һЩ����
%   % ������߼ʣ���ֹclose��ʱ��ϸ������
%   bwbord=imclearborder(bw,4);
%   % ��̬ѧ��������ʹԭ���ֿ���ϸ��ճ�ϣ���������ʹ��   
% 	se=strel('diamond',1);
%   bw_m=imclose(bwbord,se);
    bw_m = bwmorph(bw, 'majority', iteration_num);     % �����˲�����medfilt2һ��
    % ����϶
    bwfill = imfill(bw_m, 'holes');
    
%     ��������С��ǰ������һ���Ѿ��ں���������߼ʳ��ȹ��̵�ǰ�����棬ʹ���������һ���õ�ѡ��
%     STATS=regionprops(bwfill,'all');
%      for ind_i=1:length(STATS)
%          if(STATS(ind_i).Area<200)     %������С��200�ĵ�  ��Ҫ����
%              tmp=STATS(ind_i).PixelList;
%              for j=1:numel(tmp)/2
%                  bwfill(tmp(j,2),tmp(j,1))=0;
%              end
%          end
%      end

    %% 5.��ȡ��Ե
    edgeim = bwperim(bwfill); %figure;imshow(edgeim);   
    edgeim = filledgegaps(edgeim, 3);
    
    %% 6. ��Ҫ�޳�ë�̲��ܽ���edgelink�����뱣֤����Ϊ�������� (��δ�����ϸ������ʱ�����
    % ����padarray�����޸� 2015.4.28
    % �����ڲ�����
    edgeim = Remove_Spurs( edgeim );
    % ��ȥ��С������
    [LL, nn] = bwlabel(edgeim,8);
    for ii=1:nn
        if sum(sum(LL==ii))<=remove_small
            edgeim(LL==ii) = 0;
        end
    end
    % ����edgeimͼƬ
    lunkuo = edgeim; % ����������ⲿ�ռ��б���
    
    %% 7. ����������ʱ�뷽������������ʵ���ϣ�����һЩ��ֵ���״�����ӷ�ʽ�п��ܳ��ִ��ң�����˳ʱ���8�������ӣ�
    % �������edgelink�����ĸ����ԣ������޸����bug�����ѡ���ں���Ĳ����м���ֲ�
    % ����ط���10����Ŀǰ��û���õ�����Ϊ��Щbug��������Ҫ�Լ�ȥ�����̵ı߼�
    [rj, cj, ~, ~] = findendsjunctions(edgeim);
    if ~isempty(rj)
        if iteration_num==5
            for kk=1:numel(rj)
                edgeim(rj(kk),cj(kk)) = 0;
            end
        else
            lunkuo = 'error';
            e = [];
            return;
        end
    end
    [edgelist] = CX_edgelink( edgeim );
    % ��ʵmatlab�Դ����� bwtraceboundary �Ϳ���ʵ��������ܣ���2015.12.23
%     [edgelist, ~] = edgelink( edgeim, 10 );  % edgelist��ÿ��ǰ���ı�Ե���ؼ��ϣ�labelededgeim�����ص��Ӧ��ǰ�����
    
    %% 8. ȥ����С�ı߼ʣ���ο���Ϊcellɾ����ģ�壩������˫wihleѭ��ɾ��cell��ʵ�ֱ䳤��ѭ����   
    % �����ڲ����� 2016.3.9�����ɾȥ���ŵ�6.����ɾ������
%     edgelist = Remove_Small_Edgelist( edgelist, remove_small); % ȥ���߼�С�� remove_small ��ǰ��

    %% 9.��ÿ��ǰ���������зֶ��������
    seglist = lineseg( edgelist, tolerance );    %%ʹ�����Էֶ�    
    figure;
    imshow(edgeim); % �����Ļ�ͼ������edgeimǰ�������Ͻ��в���
%     imshow(zeros(size(edgeim)));
    hold on;
    
    %% ��һ���������ܽ�
%##########################################################################  

 % ����9������Ϊ��һ���̣�����ͼƬ���롢Ԥ������ȡ�����������������ֶ��������
 % ����������ڶ����̣�����Բ��ϲ���
 % ��ÿ��ǰ���������䰼����Ŀ�Ķ��٣����Ϊ 1 ��������Բ����Ҳ������ű�����Ҫʵ�ֵĹ���
 
%##########################################################################

%% ========================= ��ʼ��ÿ��ǰ����ѭ��  ========================= %%
% �������
n = length(seglist);% nΪǰ����Ŀ
% theta = cell(n, 1);    % �Ƕ���Ϣ��ʵ����û���õ���
r = {};        % ������
e = {};        % ��Բ

for i_f= 1 : n % ǰ����Ű����ϵ���ɨ��ķ�ʽ
    
%         disp('-------------------');
%         disp(['��',num2str(frame),'֡ ���ǰ��',num2str(i_f),'...']);
        % ʹ�ö��ַ���ʱ�������� seglist �� edgelist
        % ------------------ %
        ed = edgelist{i_f};
        if isequal(ed(end,:), ed(1,:)) % ͨ��ed�����һλ�͵�һλ��ͬ
        	ed = ed(1:end-1,:); % ��Ҫ�����һλȥ������ֹ�ظ�
        end
        % ------------------ %
        seg = seglist{i_f}; % seg�е����һ��Ԫ��=��һ��
        if ~isequal(seg(end,:), seg(1,:))
        	seg = [seg; seg(1,:)]; % ������ȣ����貹��һλ
        end
       %% 1.������������Ƿ���ȷ��ʹ���Լ�д��edgelink�󣬿�����֮ǰ���ж��Ƿ��д���2015.11.5��
        max_aodian = 3; % �����������ֵ���఼������Ҫ�������ݼ��ĸ��ӳ̶Ⱦ�����
        [ seg, ed, flag_aodian, hasError ] = Analyse_aodian( seg, ed );
%         if hasError
%             lunkuo = 'error';
%             close('1');
%             return;
%         end
 
       %% 4.���½�����Բ��ϲ���
        na = sum(flag_aodian);   % naΪ��������Ŀ na = ��num_aodian��
        
       %% 4.1 ����Ϊ1��0�����ȫ��
        if na<2        
            r{i_f,1} = ed;
        end
        
       %% 4.2 ����������1�����������϶��
        if na>=2
            % ע��ÿ��������Ҫ��ʼ��������ѭ��ʱ�����׳���
            % xy_axis �ǰ��������
            xy_axis=seg(flag_aodian,:);  
            % i_in_e�ǰ�����edgelist�е����
            i_in_e = zeros(na,1);      
            for jj=1:na
                i_in_e(jj) = find(ismember(ed ,xy_axis(jj,:),'rows'));
            end
            % ȥ�������İ��㣬��ѹ���ܱ�7������(��ȥ��Ҳû��ϵ 2015.10.8)
%             [ ind_in_edgelist ] = CX_sup_aodian(numel(ed)/2, ind_in_edgelist, 7);
            % ���°�����Ŀ
            na = numel(i_in_e); 

%             % ���Ƹ��º�İ���,��һ�������� ����ɾȥ��
%             x00 = ed(i_in_e,1);
%             y00 = ed(i_in_e,2);
%             plot(y00, x00, 'bs', 'LineWidth', 1.5);

            % ����ΰ�����ֶΣ�Ȼ��������ԭͼ��������
            for j2=1:na-1
                r{i_f,j2} = ed(i_in_e(j2):i_in_e(j2+1),:);    
            end
            % ǰna��Ϊ�������
            r{i_f,na} = [ed(i_in_e(na):end,:);ed(1:i_in_e(1),:)];  

            % �����ϵ����ڲ�����
            % �������cell���ɣ�����Ҫ������r����
            [ r flag_zuhe ] = Assign_PixelList( r, i_f, na ); % ��r����pixellist
            
        end % end if na>=2
        
       %% 4.3 ʹ�������������Բ
        [~, width]=size(r(i_f,:));
        % ����ÿ�е�cell��Ŀ����ǰ���������Ŀ
        while isempty(r{i_f,width}) 
                width=width-1;
        end     
            
        for k=1:width
            % ��Ϸ������Ϻã����ú��� Prasad_ElliFit
            [~, a, b, x0, y0, alpha, status] = Prasad_ElliFit(r{i_f,k}(:,2), r{i_f,k}(:,1));
            if status && a<=b*5 && size(r{i_f,k},1)>=8;    % ��ͼ��Լ��������������̫��ֵ���Բ
                c=cos(alpha);
                s=sin(alpha);
                polar_angle=linspace(0,360,361);
                xq= a*cosd(polar_angle);
                yq= b*sind(polar_angle);
                xn=xq*c-yq*s+x0;
                yn=xq*s+yq*c+y0;
                plot(xn,yn,'r','LineWidth',1.5);
                % ������Բ����,������Բʲô���������棨���ڿ��޸�Ϊֻ����status��Ϊ0�ģ�
                e{i_f,k}.a=a;
                e{i_f,k}.b=b;
                e{i_f,k}.x0=x0;
                e{i_f,k}.y0=y0;
                e{i_f,k}.alpha=alpha/pi*180; % ͼƬ��y�������µģ���˽Ƕ���ʱ��Ϊ-������������ϵ����ʱ��Ϊ+
                e{i_f,k}.status=status;
                e{i_f,k}.similar_delete = 0; % ���������ԣ�Ϊ1����ɾ����������
                e{i_f,k}.hd = 0;
                e{i_f,k}.num_hypoth = width;
                e{i_f,k}.ind_region = i_f;
                % ������������
                if width>1
                    e{i_f,k}.flag_combine = flag_zuhe(k,:);
                else
                    e{i_f,k}.flag_combine = 1;
                end
            else % ������������statusΪ0
                e{i_f,k}.status = 0;
            end

        end % һ��ǰ���������
        
       %% 4.4 �����˹������
%        disp(['������˵',num2str(width),'��']);
%        disp('Calculate Hausdorff distance...')
       e(i_f,:) = CX_Hausdorff_Distance(e(i_f,:), r(i_f,:));  % �����˹�������
%        disp('��ɣ�');

end
hold off;

toc;

end



% ==================== ���²���Ϊ CX_FIT �ڲ����õĺ��� =================== %
%
% ��Щ����ԭ�����ⲿ���������룬�������
% ȫ������path֮���Ƿ�������Ƚ����� 2015.10.7
% ======================================================================= %





















