function [im, count] = plot_ellipse_label(im, Label, e, scale)

% % ˵���ǵ���˵ǰ����������ǰ����״��������Բ label����Ϊ0�������Ϊ������
ny0 = round(e.y0); ny0(ny0<=0) = 1; ny0(ny0>size(im,1)) = size(im,1);
nx0 = round(e.x0); nx0(nx0<=0) = 1; nx0(nx0>size(im,2)) = size(im,2);
if e.num_hypoth==1 && Label(ny0, nx0)~= 0 
    logi = Label==Label(ny0, nx0);
    count = sum(sum(logi));
    % ------------------------------------------------------------------- %
    % ����ط���bug���ڣ���2��ʮ�ֽ����������Բ���ʱ�Ƕ����ģ���bwlabelʱ����ͬһ��
    % ʹ��һ����Բռ��2�����򣬵��³���
    % The automatic track with label 49 is not consistent with the image data!
    % �����Ҫ���������ֵ�жϣ���������ݼ�������seglistʱ�������ִ���
    % ����������������Բ���̫�࣬�򲻲��������ʾ
    % ------------------------------------------------------------------- %
    if count && count<=1.5*pi*e.a*e.b % count����Ϊ0
        im(logi) = e.color;
        return % �����������ڼ�ǿ�ָ��
    end
end

count = 0;
e.a = e.a*scale; % ��΢��С����Բ
e.b = e.b*scale;


c=sqrt(e.a^2-e.b^2);   % ����
cc=cosd(e.alpha);
ss=sind(e.alpha);
c1=[c*cc+e.x0,c*ss+e.y0];
c2=[-c*cc+e.x0,-c*ss+e.y0];

squares = 5;
for x = floor(e.x0-e.a) : ceil(e.x0+e.a)
    if x<1 || x>size(im,2)
        continue
    end
    for y = floor(e.y0-e.a) : ceil(e.y0+e.a)
        if y<1 || y>size(im,1) % ����Ĳ���
            continue
        end
        % Ҫ��õ���벻�ܱ����labelռ��
        if norm([x,y]-c1) + norm([x,y]-c2) <= 2*e.a && im(y,x)==0
            % ������Բ�ڵĵ�Ҷ�ֵ
            im(y,x) = e.color;
            count = count + 1;
        end
        if count==0 && y==round(e.y0) && x==round(e.x0) % ����Ҳ�����ŵ㣬�򸲸�һ������ȷ���ܼ���AOG��TRA
            im(y-squares:y+squares,x-squares:x+squares) = e.color*ones(squares*2+1);
            count = count + (squares*2+1)^2;
        end
    end
end






