function  e = CX_Ellipse_Optimal( e )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% ����һ�п�����debug
% e = Ellipse;
frame = numel(e);

for t=1:frame 
    disp(['�Ż�ɸѡ��',num2str(t),'֡...']);
    flag = ~isemptycell(e{t});
    height = size(flag,1);
    for hh=1:height        %%ÿ��ǰ��
        
        width = sum(flag(hh,:));  %%����ÿ�е�cell��Ŀ����ǰ���������Ŀ
        
        %% ��һ����ɾȥstatusΪ0�ļ�˵
        %########### ɾȥstatusΪ0�ļ�˵ ##########
		tmp_e = e{t}(hh,:);
		for i=1:width
			if tmp_e{i}.status==0
				tmp_e{i} = [];
			end
		end
		sts_flag = ~isemptycell(tmp_e);
        new_width = sum(sts_flag); % ɾ����ĳ���
        if new_width==0 % ���ж�ɾû��
            e{t}(hh,1) = {[]};
            continue
        end  
        e{t}(hh,1:new_width) = tmp_e(sts_flag);
        e{t}(hh,new_width+1:end) = {[]};   

        %% �ڶ������ڶ఼���������Բ��˵̫�࣬ɾ��һЩ��ϳ̶Ƚϲ��
         %######## ��ֺ���������Ҫ����######## 
        tmp_e = tmp_e(sts_flag);
        % B. ֱ����cell������������,����
        tmp_fitness = zeros(1,new_width);
        for i=1:new_width
            % ��ֺ����ݶ�Ϊ hd/sqrt(num_pixels)
            tmp_fitness(i) = tmp_e{i}.hd/sqrt(tmp_e{i}.num_pixels);
        end
        [~, ind_fit] = sort(tmp_fitness);
        sort_e = tmp_e;
        for i=1:new_width
            sort_e(i) = tmp_e(ind_fit(i));
        end
        % �������һ��ǰ�����Խ���ɾ�����ٷ���ԭcell
        if new_width >15
            num_remain = 15; % ��˵����10��������϶���õ�10��������ȫ������
        else
            num_remain = new_width; % ���屣���ĸ���
        end
        % ��e��ֵ
        e{t}(hh, 1:num_remain) = sort_e(1:num_remain);
        e{t}(hh, num_remain+1:end) = {[]};
        
        %% ��������Ѱ��������Բ ��2015.6.27�޸���һ��bug��������� new_width ����Ϊ num_remain ��
        for ind_1=1:num_remain-1
            for ind_2=ind_1+1:num_remain
                e1=e{t}{hh,ind_1};
                e2=e{t}{hh,ind_2};
                % �����Բ�޽�����ֱ�ӽ�����һ��ѭ��
                if max([abs(e1.x0-e2.x0), abs(e1.y0-e2.y0)]) > e1.a + e2.a
                    continue;
                end

%                 e1_4_c = struct2array(e1);
%                 e2_4_c = struct2array(e2);   
                % �·�����������е����⣬�޷���ȷ����������Ȳ�����
%                 overlap = CX_calculate_overlap_C(e1_4_c,e2_4_c);%%ʹ��C��д���ص��Ⱦ���72�룬��Ϊ������������17��

                overlap = CX_calculate_overlap(e1,e2);   % CX_calculate_overlap ��matlab��ļ�����̣�166��
                if overlap >= 0.75    %%�ص��ȶ���Ϊ0.9
%                     belta_i = e1.num_pixels;
%                     belta_j = e2.num_pixels;   %%ֱ�ӱȽ���ϱߵĳ���
                    e{t}{hh,ind_2}.similar_delete = 1; % ��Ϊe2���ں��棬˵����ϳ̶Ƚϲ��˽���ɾ��
                    
                if 0 % ע�͵��ķ���
                  %####  �ж�Բ�Ľǵķ����ٶ�̫��
    %             %%�ж�Ϊͬһ��ʱ��Ƚ�Բ�ĽǴ�С
    %                 rr=r{n,j}*[0 1;1 0];   %%rr(y,x)��Ҫ���в��ܱ��rr(x,y)
    %                 v1=rr(1,:)-[e2.x0,e2.y0];
    %                 v2=rr(round(end/2),:)-[e2.x0,e2.y0];
    %                 v3=rr(end,:)-[e2.x0,e2.y0];
    %                 theta1=acosd(v1*v2'/sqrt(v1*v1'*v2*v2'));
    %                 theta2=acosd(v2*v3'/sqrt(v2*v2'*v3*v3'));
    %                 belta_j=theta1+theta2;
    %                 temp_e{n,i}.theta=belta_i; 
    
%                     if belta_i>=belta_j
%                         e{t}{n,j}.similar_delete=1;   %%����ϱ߳���С���Ǹ�similar_delete=1���ں����Ͳ�������
%                     else
%                         e{t}{n,i}.similar_delete=1;
%                     end
                end
                
                end
            end
        end
        
        %% ���Ĳ���ɾ��������Բ������ num_hypoth �����жϸ�ǰ���Ƿ�Ϊ����˵ǰ��
        % --------------------------------------------------------------- %
        % �����ⲿ�ֺ���Ҫ�������� num_hypoth �Ļ��ᵼ�º����˹����ʱ����
        % ������Ҫ���� num_hypoth �жϸ���Բ�Ƿ�λ�ڶ��˵ǰ����
        % --------------------------------------------------------------- %
        % ���ɾȥ similar_delete=1 ����Բ��������ͳ�� num_hypoth  2015.6.28
		tmp_e = e{t}(hh,:);
		for i=1:num_remain
			if tmp_e{i}.similar_delete==1
				tmp_e{i} = [];
			end
		end
		simd_flag = ~isemptycell(tmp_e);
        num_remain = sum(simd_flag); % ɾ����ĳ���
		e{t}(hh,1:num_remain) = tmp_e(simd_flag);
        e{t}(hh,num_remain+1:end) = {[]};

        % �����ǰ����������Բ�� num_hypoth
        for ii=1:num_remain
            e{t}{hh,ii}.num_hypoth = num_remain;
        end
        
        %% ���岽������ sum��flag_combine���������򣬽�������ϵķ���ǰ��
        % --------------------------------------------------------------- %
        % �ܽ� 2015.7.7 ---
        
        % ���е�����ʱ��ÿ��ǰ���е� status Ϊ0�� similar_delete Ϊ1����Բ�Ѿ�
        % ɾȥ�����Ұ�����ϳ̶�ȡǰ10λ�������򣨴���10���Ľ�ȡ10����С��10����ȡԭ��Ŀ��
        %
        % Ȼ���ں�������ì�ܼ�˵ʱ����Ϊ�Ǵ�ǰ����ȡ numel(flag_combine) ����1��Բ��Ϊ����
        % ���������е� 100 010 001������չ��3��ì�ܼ���
        % 
        % ��ˣ� ������Ĵ�������Ҫ��ÿһ�е���Բ������������ʹ�õ�����ǰ�档
        % *******��ע�⣬���������ûʵ��ʱʹ�õ��� 22ì�ܵļ�˵���������Լ���������࣬Ӱ������ٶȣ� ******
        
        tmp_e = e{t}(hh, 1:num_remain);
        sum_flag_combine = zeros(1, num_remain);
        % ��� flag_combine �ĺ�
        for ind_e=1:num_remain
            sum_flag_combine(ind_e) = sum( tmp_e{ind_e}.flag_combine );
        end
        % ����С�ķ���ǰ�棬�����η�ǰ��
        [ ~, new_ind ] = sort(sum_flag_combine);
        % -------- %
        % ��ֵ����e�����̽���
        e{t}(hh, 1:num_remain) = tmp_e(new_ind);
        e{t}(hh, num_remain+1:end) = {[]};
        % --------------------------------------------------------------- %
        
    end
    %% �����Ż����ԣ����ޣ�
      
end

end % end function


%% �Ѿ���дΪC�汾 �������matlab�汾�����ˣ���������c�汾�����⣬�����ʹ�����棩

function [ overlap ] = CX_calculate_overlap(e1,e2)

        c1=sqrt(e1.a^2-e1.b^2);   %%����
        cc1=cosd(e1.alpha);
        ss1=sind(e1.alpha);
        c11=[c1*cc1+e1.x0,c1*ss1+e1.y0];
        c12=[-c1*cc1+e1.x0,-c1*ss1+e1.y0];
        
        count1 = pi*e1.a*e1.b;
        count2 = pi*e2.a*e2.b;
        count=0;
        
        c2=sqrt(e2.a^2-e2.b^2);   %%����
        cc2=cosd(e2.alpha);
        ss2=sind(e2.alpha);
        c21=[c2*cc2+e2.x0,c2*ss2+e2.y0];
        c22=[-c2*cc2+e2.x0,-c2*ss2+e2.y0];

        for x=e2.x0-e2.a : e2.x0+e2.a
            for y=e2.y0-e2.a : e2.y0+e2.a
                if sqrt((x-c21(1))^2+(y-c21(2))^2)+sqrt((x-c22(1))^2+(y-c22(2))^2)<=2*e2.a
                    if sqrt((x-c11(1))^2+(y-c11(2))^2)+sqrt((x-c12(1))^2+(y-c12(2))^2)<=2*e1.a
                        count=count+1;
                    end
                end
            end
        end
        overlap=count/(count1+count2-count);
end

