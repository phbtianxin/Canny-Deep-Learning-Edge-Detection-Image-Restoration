function smoV = testSmooth_f(image)
%��һ��2ֵ����ͼ���е��߶���ϳ�һ�����ߣ���ͨ��������ϲ������ʵ�ƽ���仯�ʵõ������۳̶�

%����ʾ����
% smoV = testSmooth_f('line04.png');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%��ȡ2ֵ����ͼ��
ima = imread(image); 
ima = ima(:,:,1);

%��ʼ��
[r,c] = find(ima==0);
xy(:,1) = c; %��������Ϊ������
xy(:,2) = r; %��������Ϊ������

%��ʱ�Զ�Ѱ����϶�
t = [xy(2:end,1);xy(end,1)];
d = t - xy(:,1);
gap = find(d>3);
p1 = [xy(gap,1),xy(gap,2)]; %��������ӵ�һ��
p2 = [xy(gap+1,1),xy(gap+1,2)]; %��������ӵ���һ��

%���������������������ߵ���������


%�Ƕȹ�һ��
v = [p2(1)-p1(1),p2(2)-p1(2)];
angle = atan(v(2)/v(1));
newXY = [xy(:,1),xy(:,2)]*[cos(angle),-sin(angle);sin(angle),cos(angle)];
newp1 = p1*[cos(angle),-sin(angle);sin(angle),cos(angle)];
newp2 = p2*[cos(angle),-sin(angle);sin(angle),cos(angle)];

%��ֱλ�ù�һ��
newp1(:,2) = newp1(:,2) - newXY(gap,2);
newp2(:,2) = newp2(:,2) - newXY(gap,2);
newXY(:,2) = newXY(:,2) - newXY(gap,2);

%ˮƽλ�ù�һ��
newp1(:,1) = newp1(:,1) - min(newXY(:,1));
newp2(:,1) = newp2(:,1) - min(newXY(:,1));
newXY(:,1) = newXY(:,1) - min(newXY(:,1));
newXY(1:find(newXY(:,1)==0)-1,:) = []; %��ȥ��˵��ظ�Ӱ�䲿��
newXY(find(newXY(:,1)==max(newXY(:,1)))+1:end,:) = []; %��ȥ��˵��ظ�Ӱ�䲿��

%%���ƹ�һ�������ɢ������
figure(1)
plot(newXY(:,1),newXY(:,2))
title('Unitary discrete curve','FontSize',10,'FontAngle','italic');

%���
f = fit(newXY(:,1),newXY(:,2),'poly5');

%��������ϲ������ʵ�ƽ���仯�������߹⻬��
cc(1) = f.p1;
cc(2) = f.p2;
cc(3) = f.p3;
cc(4) = f.p4;
cc(5) = f.p5;
cc(6) = f.p6;
x=sym('x');
cf = cc*[x^5;x^4;x^3;x^2;x;1]; %���ߺ���
cK = abs(diff(diff(cf)))/(1+(diff(cf))^2)^(3/2); %���ߵ����ʺ���
dcK = diff(cK); %�������ʵĵ�����
dx = floor(newp1(1,1))+1:floor(newp2(1,1)); %������
SdcK = abs(subs(dcK,x,dx)); %�������㴦�����ʱ仯��
smoV = sum(SdcK)/length(dx); %������ϲ��ֵ�ƽ�����ʱ仯��

%%�����������
X = min(newXY(:,1)):max(newXY(:,1));
Y = (f(X))';
figure(2),
plot(newXY(:,1),newXY(:,2),'*')
hold on
plot(X,Y,'r')
title('Fitted curve','FontSize',10,'FontAngle','italic');

%%������ϲ��ֵ����ʺ���ͼ��
ScK = subs(cK,x,dx);
figure(3),
plot(dx,ScK)
title('Curvature function','FontSize',10,'FontAngle','italic');

%%������ϲ��ֵ����ʱ仯�ʾ���ֵ����ͼ��
figure(4),
plot(dx,SdcK)
title('Curvature absolute derivative function','FontSize',10,'FontAngle','italic');

