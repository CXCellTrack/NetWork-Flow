function [ e_row ] = CX_Hausdorff_Distance( e_row, r_row )
%% ������Ŀ����Բ������ϵ��
% ������Ŀ    ��Բ����
%   1           1
%   2           3
%   3           7
%   4           15
%   5           31
%   k          2^k-1

%%

flag = ~isemptycell(e_row);
% for n=1:size(e,1) % ÿ��ǰ��
width = sum(flag); % ǰ���еļ�˵��Ŀ
    
for i=1:width
    if e_row{i}.status==0 % status Ϊ0 ������
        continue;
    end
    e_row{i} = calculate_distance(e_row{i} ,r_row{i}); % ���뵱ǰԲ����ǰ�㼯����ǰȫ��r{n,width}���ƺ��ò��ϣ�
end

% end
end



function  [ e_n_i ]= calculate_distance( e_n_i ,r_n_i )  %%����Ҫ��ֵ���ݳ��������ܵ�ַ����

    [ exy ] = ellipse_xy(e_n_i);
    rr = r_n_i*[0 1;1 0];        %%����������2�н���   r�Ǽ̳���edgelist,seglist��2������xy�෴�ģ���Ҫ��λ
%     r_all = r_n_width*[0 1;1 0];   %%r_all��ʾ��������
    [ hd ~ ] = HausdorffDist(rr, exy, 1);    %%�˺����Ѿ�����Ϊ������Բ��Χ�����������̾��룬�����϶ȸ���
    
    %% ����g��cl��e����Ŀǰ�е����⣬�Ȳ���
%     complement_rr=setdiff(r_all, rr, 'rows');  %%rr�Ĳ��������ڼ���g��cl��e��
%     flag= false(numel(complement_rr)/2,1);  %%��־����
%     for j=1:numel(complement_rr)/2
%         xg=complement_rr(j,1);
%         yg=complement_rr(j,2);
%         if sqrt((complement_rr(j)-c11)*(complement_rr(j)-c11)')+sqrt((complement_rr(j)-c12)*(complement_rr(j)-c12)')<=2*ee.a  %%��ȡrr�Ĳ�����λ����Բ�ڵĵ�
%             flag(j)=1;
%         end
%     end
%     CL=complement_rr(flag,:);
%     if ~isempty(CL)
%         [gd ~] = HausdorffDist(CL,exy,1,'vis');   %%����g��cl��e��,�վ���Ͳ����ˣ��ӿ��ٶ�
%     else
%         gd = 0;
%     end
    %%
%     gd=0;
    e_n_i.hd = hd;
    e_n_i.num_pixels = size(rr,2); % 2015.5.25
%     e_n_i.shape_yueshu = shape_yueshu;
%     e_n_i.Cl_distance = (hd + gd)*shape_yueshu;  %%���빫ʽ
    
end

function [ exy ]= ellipse_xy(e_n_i)   %%������Բe{i,j}�������Բ�ϵĵ㼯�͹�ʽ��2��,����2������
    % ������Բ���ϵĵ�����
    c=cosd(e_n_i.alpha);
    s=sind(e_n_i.alpha);
    polar_angle=linspace(0,360,181);
    xq= e_n_i.a*cosd(polar_angle);
    yq= e_n_i.b*sind(polar_angle);
    xn=xq*c-yq*s+e_n_i.x0;
    yn=xq*s+yq*c+e_n_i.y0;
    exy=[xn',yn'];    %%�õ���Բ���ϵĵ�

%     Ce=pi*sqrt(2*(ee.a^2+ee.b^2));    %%��Բ�ܳ����ƹ�ʽ
%     ec=sqrt(ee.a^2-ee.b^2)/ee.a;   %%ƫ��   c/a
%     shape_yueshu=Ce*sqrt(1/(1+ec^2));
%     
% 	ccc=sqrt(ee.a^2-ee.b^2);   %%�뽹��c
%     c11=[ccc*c+ee.x0 ,ccc*s+ee.y0];
%     c12=[-ccc*c+ee.x0 ,-ccc*s+ee.y0];

    
end
 

    
    
    
