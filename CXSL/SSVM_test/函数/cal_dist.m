function method = cal_dist(s1,s2,e1,e2)

% 进行匹配，找到相近的作为迁移
d11 = sqrt((s1.x0-e1.x0)^2 + (s1.y0-e1.y0)^2);
d22 = sqrt((s2.x0-e2.x0)^2 + (s2.y0-e2.y0)^2);

d12 = sqrt((s1.x0-e2.x0)^2 + (s1.y0-e2.y0)^2);
d21 = sqrt((s2.x0-e1.x0)^2 + (s2.y0-e1.y0)^2);

if d11+d22<=d12+d21
    method = 1; % d=1表示正常匹配，1对1，2对2
else
    method = 0; % 否则1对2，2对1
end