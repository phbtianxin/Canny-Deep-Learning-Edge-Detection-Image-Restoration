function smoV = testSmooth_f(image)
%将一幅2值曲线图像中的线段拟合成一条曲线，并通过求其拟合部分曲率的平均变化率得到其曲折程度

%调用示例：
% smoV = testSmooth_f('line04.png');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%读取2值曲线图像
ima = imread(image); 
ima = ima(:,:,1);

%初始化
[r,c] = find(ima==0);
xy(:,1) = c; %以列数作为横坐标
xy(:,2) = r; %以行数作为纵坐标

%暂时自动寻找拟合端
t = [xy(2:end,1);xy(end,1)];
d = t - xy(:,1);
gap = find(d>3);
p1 = [xy(gap,1),xy(gap,2)]; %需拟合连接的一端
p2 = [xy(gap+1,1),xy(gap+1,2)]; %需拟合连接的另一端

%按曲线走向重新排序曲线的坐标序列


%角度归一化
v = [p2(1)-p1(1),p2(2)-p1(2)];
angle = atan(v(2)/v(1));
newXY = [xy(:,1),xy(:,2)]*[cos(angle),-sin(angle);sin(angle),cos(angle)];
newp1 = p1*[cos(angle),-sin(angle);sin(angle),cos(angle)];
newp2 = p2*[cos(angle),-sin(angle);sin(angle),cos(angle)];

%垂直位置归一化
newp1(:,2) = newp1(:,2) - newXY(gap,2);
newp2(:,2) = newp2(:,2) - newXY(gap,2);
newXY(:,2) = newXY(:,2) - newXY(gap,2);

%水平位置归一化
newp1(:,1) = newp1(:,1) - min(newXY(:,1));
newp2(:,1) = newp2(:,1) - min(newXY(:,1));
newXY(:,1) = newXY(:,1) - min(newXY(:,1));
newXY(1:find(newXY(:,1)==0)-1,:) = []; %砍去左端的重复影射部分
newXY(find(newXY(:,1)==max(newXY(:,1)))+1:end,:) = []; %砍去左端的重复影射部分

%%绘制归一化后的离散点曲线
figure(1)
plot(newXY(:,1),newXY(:,2))
title('Unitary discrete curve','FontSize',10,'FontAngle','italic');

%拟合
f = fit(newXY(:,1),newXY(:,2),'poly5');

%根据其拟合部分曲率的平均变化率求曲线光滑度
cc(1) = f.p1;
cc(2) = f.p2;
cc(3) = f.p3;
cc(4) = f.p4;
cc(5) = f.p5;
cc(6) = f.p6;
x=sym('x');
cf = cc*[x^5;x^4;x^3;x^2;x;1]; %曲线函数
cK = abs(diff(diff(cf)))/(1+(diff(cf))^2)^(3/2); %曲线的曲率函数
dcK = diff(cK); %曲线曲率的导函数
dx = floor(newp1(1,1))+1:floor(newp2(1,1)); %采样点
SdcK = abs(subs(dcK,x,dx)); %各采样点处的曲率变化率
smoV = sum(SdcK)/length(dx); %整个拟合部分的平均曲率变化率

%%绘制拟合曲线
X = min(newXY(:,1)):max(newXY(:,1));
Y = (f(X))';
figure(2),
plot(newXY(:,1),newXY(:,2),'*')
hold on
plot(X,Y,'r')
title('Fitted curve','FontSize',10,'FontAngle','italic');

%%绘制拟合部分的曲率函数图像
ScK = subs(cK,x,dx);
figure(3),
plot(dx,ScK)
title('Curvature function','FontSize',10,'FontAngle','italic');

%%绘制拟合部分的曲率变化率绝对值函数图像
figure(4),
plot(dx,SdcK)
title('Curvature absolute derivative function','FontSize',10,'FontAngle','italic');

