%最后用我最熟悉的matlab给你画个圣诞树吧
clc

clear

close all

figure('color',[0 0 0])

title('Merry Christmas','color','w')
hold on
t = 10*pi:-0.1:0;
A = 1;w = 1;
for i = 1:length(t)   
    x(i) = A *cos(w*t(i));   
    y(i) = A *sin(w*t(i));    
    A = A + 0.1;
end
z = t*20;

scatter3(x,y,z,'g.','MarkerFaceColor',[0 1 0]);

hold on
%红色螺旋曲线

scatter3(-x,-y,z,'r.','MarkerFaceColor',[1 0 0]);

axis off

%随机生成200个雪花
snow=200;snowX=(-50 + (50+50) * rand(snow,1));
snowY=(-50 + (50+50) * rand(snow,1));
snowZ=(rand(snow,1)*800);
scatter3(snowX,snowY,snowZ,'w*')
%圣诞树顶的黄色五角星

p = plot3(0,0,z(1)+20,'Pentagram','MarkerSize',20,'MarkerEdgeColor','y','MarkerFaceColor','y');
hold off



%旋转50次，不同的角度观察圣诞树

for i = 1:50

view(-1*i,0)

i = i+1;

pause(0.1)

end


hold off

