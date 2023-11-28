function BObj = optimizeByGA_f(edgPm,M,N)
%���ݲ�����Եĸ������ߵĲ���ͨ��GA�㷨�������ƥ�䷽��

%ǰ����裺
% 1��һ��������һ�����������
% 2����Ϊֻ�д������������������ڳɻ�����Ƶ��������߲���ƥ���ʸ�

%���У�
% BObj -- һ��m+1ά��������m��������������������ӵ�һ��Ԫ�ؿ�ʼ��ÿ������Ԫ��������Ҫ������Ե���������ţ����һ��Ԫ������ԵĶ�������λ��0��䡣
% edgPm -- һ��Lind��N���󣬱�ʾLind���μ����������ÿ���ߵ�N����һ��ƥ�������
% M,N -- ͼ��ĳߴ硣

%���磺
% [edgPm,w,edgs,M,N,Q] = getEdge_f('planet01.png');
% BObj = optimizeByGA_f(edgPm,M,N);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%������ʼ��%%%
Nind = 20; %������Ⱥ������
MaxGen = 100; %��������Ŵ�����
T = 1.2; %����ƥ������
Lind = size(edgPm,2); %�����ɸ���Ĳ������������μ�ƥ�������������
%BaseV = [Lind:-1:2,floor(Lind/2)]; %��������������
BaseV = [Lind:-1:1,floor(Lind/2)]; %��������������
Chrom = crtbp(Nind,BaseV); %��ó�ʼ��Ⱥ
Pm = 0.7/Lind; %���ñ������
XOVR = 0.3; %���ý������

traceM = []; %��ʼ��һ����Ⱥƽ�����ܸ�������
traceB = []; %��ʼ��һ����Ⱥ������ܸ�������

%%%�����ÿ���������ߵ�ƥ�����%%%
difference = zeros(Lind); %��ʼ��һ��Lind�������ߵ��������׷���������¼ÿ����������֮���ƥ�����
    
%��������������۳̶�
figurek = 1; %��ʼ����ͼ���
Tdifference = 0*difference; %��ʼ��һ����ʱ��¼����
for i=1:Lind-1
    for j=i+1:Lind %ֻ���õ�difference������ϰ�Ǽ���,���������Խ���
        
        %��ʼ��
        xy1 = [edgPm(i).Curve(:,2),edgPm(i).Curve(:,1)]; %�ɵڶ˵�1�����߶ε������к������зֱ�õ�����1��x�������к�y��������
        xy2 = [edgPm(j).Curve(:,2),edgPm(j).Curve(:,1)]; %�ɵڶ˵�2�����߶ε������к������зֱ�õ�����2��x�������к�y��������
        p1 = [edgPm(i).Beginning(2),edgPm(i).Beginning(1)]; %�õ�һ�����𴦶˵㣨�˵�1���ĺ�������
        p2 = [edgPm(j).Beginning(2),edgPm(j).Beginning(1)]; %�õ���һ�����𴦶˵㣨�˵�2���ĺ�������
        
        %%������ɢ������
        %figure(figurek)
        %subplot(2,2,1)
        %plot(xy1(:,1),xy1(:,2),'.',xy2(:,1),xy2(:,2),'.')
        %title('Original discrete curve','FontSize',10,'FontAngle','italic');
        %axis equal %ʹͼ��ĺ����������ʾ����һ��
        
        %�Ƕȹ�һ��
        v = [p2(1)-p1(1),p2(2)-p1(2)]; %���˵������
        angle = atan(v(2)/v(1)); %���˵�����������
        newXY1 = xy1*[cos(angle),-sin(angle);sin(angle),cos(angle)]; %�����������ת����1����������
        newXY2 = xy2*[cos(angle),-sin(angle);sin(angle),cos(angle)]; %�����������ת����2����������
        newp1 = newXY1(1,:); %�������ת�Ķ˵�1����
        newp2 = newXY2(1,:); %�������ת�Ķ˵�2����
        
        %��������1���ظ�ӳ�䲿��
        mean5XY1 = (newXY1(5:end,1)+newXY1(4:end-1,1)+newXY1(3:end-2,1)+newXY1(2:end-3,1)+newXY1(1:end-4,1))/5; %�õ���һ������x�����5�����
        dM5XY1 = mean5XY1(1:end-1) - mean5XY1(2:end); %����֡��õ�ÿ��������ֵ����һ��������ֵ�Ĳ�ֵ
        %position1 = find(dM5XY1(1)*dM5XY1<0); %�����������ֵ���Ÿı��λ��
        position1 = find(round(dM5XY1)==0); %�����������ֵ��С��λ��
        if length(position1)~=0
            newXY1 = newXY1(1:position1(1)+4,:); %ֻ�����Ӷ˵㵽��һ��������5���ֵ�仯��С�ĵ�
        end
        
        %��������2���ظ�ӳ�䲿��
        mean5XY2 = (newXY2(5:end,1)+newXY2(4:end-1,1)+newXY2(3:end-2,1)+newXY2(2:end-3,1)+newXY2(1:end-4,1))/5; %�õ���һ������x�����5�����
        dM5XY2 = mean5XY2(1:end-1) - mean5XY2(2:end); %����֡��õ�ÿ��������ֵ����һ��������ֵ�Ĳ�ֵ
        %position2 = find(dM5XY2(1)*dM5XY2<0); %�����������ֵ���Ÿı��λ��
        position2 = find(round(dM5XY2)==0); %�����������ֵ��С��λ��
        if length(position2)~=0
            newXY2 = newXY2(1:position2(1)+4,:); %ֻ�����Ӷ˵㵽��һ��������5���ֵ�仯��С�ĵ�
        end
        
        %%���ƽǶȹ�һ����ȥ��Ӱ�����ɢ������
        %figure(figurek)
        %subplot(2,2,2)
        %plot(newXY1(:,1),newXY1(:,2),'.',newXY2(:,1),newXY2(:,2),'.')
        %title('Unitary discrete curve','FontSize',10,'FontAngle','italic');
        %axis equal %ʹͼ��ĺ����������ʾ����һ��
        
        %%%�������������ϵ����������
        if ~(((newXY1(end,1)<newXY1(1,1))&&(newXY1(1,1)<newXY2(1,1))&&(newXY2(1,1)<newXY2(end,1)))||...
                ((newXY1(end,1)>newXY1(1,1))&&(newXY1(1,1)>newXY2(1,1))&&(newXY2(1,1)>newXY2(end,1)))) %��������ߴ����������ⲻ��ɢ
            Tdifference(i,j) = inf; %��������߲��ʺ�ƥ�䣬�Ӷ�ƥ�����Ϊinf
            continue
        end
        
        %���
        %%%�����Զ�ѡ��������ߵĴ���
        f = fit([newXY1(:,1);newXY2(:,1)],[newXY1(:,2);newXY2(:,2)],'poly5');

        %��������ϲ������ʵ�ƽ���仯�������߹⻬��
        clear cc;
        cc(1) = f.p1; %��1��ϵ��
        cc(2) = f.p2; %��2��ϵ��
        cc(3) = f.p3; %��3��ϵ��
        cc(4) = f.p4; %��4��ϵ��
        cc(5) = f.p5; %��5��ϵ��
        cc(6) = f.p6; %��6��ϵ��
%         cc(7) = f.p7; %��7��ϵ��
%         cc(8) = f.p8; %��8��ϵ��
%         cc(9) = f.p9; %��9��ϵ��
%         cc(10) = f.p10; %��10��ϵ��
        x=sym('x'); %��x����Ϊһ������
        cf = cc*[x^5;x^4;x^3;x^2;x;1]; %������ߵķ��ź���
        cK = abs(diff(diff(cf)))/(1+(diff(cf))^2)^(3/2); %���ߵ����ʷ��ź���
        dcK = diff(cK); %�������ʵĵ����������ţ�
        dx = round(newp1(1)):(newp2(1)-newp1(1))/abs(newp2(1)-newp1(1)):round(newp2(1)); %������
        SdcK = abs(subs(dcK,x,dx)); %�������㴦�����ʱ仯��
        difference(i,j) = sum(SdcK)/length(dx); %������ϲ��ֵ�ƽ�����ʱ仯��
        
        %%������ϲ��ֵ�����
        %figure(figurek)
        %subplot(2,2,3)
        %plot(newXY1(:,1),newXY1(:,2),'.',newXY2(:,1),newXY2(:,2),'.')
        %hold on
        %plot(dx,(f(dx))','r.')
        %title('Fitted curve','FontSize',10,'FontAngle','italic');
        %axis equal %ʹͼ��ĺ����������ʾ����һ��
        
        %%������ϲ��ֵ����ʱ仯�ʺ���ͼ��
        %ScK = subs(cK,x,dx);
        %figure(figurek)
        %subplot(2,2,4)
        %plot(dx,SdcK)
        %title('The derivative of curvature function','FontSize',10,'FontAngle','italic');
        %axis equal %ʹͼ��ĺ����������ʾ����һ��
        
        figurek = figurek + 1;
    end
end
difference = difference./(sum(sum(difference))/(sum(sum(difference>0)))); %��һ���������Ը�������ƽ��ֵ��
difference = difference + Tdifference;
    
%����ƽ��ɫֵ����̶�
Tdifference = 0*difference; %��ʼ��һ����ʱ��¼����
for i=1:Lind
    for j=i+1:Lind %ֻ���õ�difference������ϰ�Ǽ���
        Tdifference(i,j) = abs(edgPm(i).MeanColor(1)-edgPm(j).MeanColor(2))+abs(edgPm(i).MeanColor(2)-edgPm(j).MeanColor(1));
    end
end
Tdifference = Tdifference./(sum(sum(Tdifference))/(sum(sum(Tdifference>0))));%��һ���������Ը�������ƽ��ֵ��
difference = difference + Tdifference; %���ƥ��ĸ��ֲ����ܺ�
difference = difference + difference'; %�����Ǿ�����������

%%%�����Ŵ��Ż�%%%
gen = 0; %���ó�ʼ����
while gen < MaxGen
    
%--------����--------%
    Chrom = mut(Chrom,Pm,BaseV);
%--------------------%
    
%--------����--------%
    Chrom = xovsp(Chrom, XOVR);
%--------------------%
    
%--------����ɱ�����--------%
    Chrom = Chrom+1; %����������1��Ϊ���
    Phenotype = 0*Chrom; %��������;��󣬳����һ���������������ÿ����0Ԫ�ر�ʾ�μ���Ե���������š��ӵ�һ��Ԫ�أ�ÿ������Ԫ����ԡ�
    for i=1:Nind
        tem = 1:Lind; %һ�����ڷ�������ŵ��б�
        for j=1:2*Chrom(i,end) %ֻ����Ҫ�μ���Ե��߱�ż���
            Phenotype(i,j) = tem(Chrom(i,j)); %�˴������ߵ�������ִ�����б��е�Chrom(i,j)λ���ϵ�ֵ
            tem(Chrom(i,j)) = []; %�ѱ�ѡ�ߵ������Ҫ��������б����Ƴ�
        end
    end
    %Phenotype = [Phenotype,Chrom(:,Lind)]; %����ÿ�ַ�������Ե������߶���
    %avoid multi column zero
    Phenotype(:,end)=Chrom(:,Lind+1);
    Chrom = Chrom-1; %�������ͻ�ԭΪ��0��Ϊ���
%--------------------%
    
%--------����ÿ�ַ�����ƥ����ۣ�Ŀ�꺯��ֵ��--------%
    tPhenotype = Phenotype + (Phenotype==0); %%
    %�������;������е�0Ԫ�ض����1���Ӷ��Ѳ���0����0���ߣ�������ԵĿ�λ����ƥ������ɲ���1����1�ŵ�ƥ����죨��õ�0��
    ObjV = zeros(Nind,1); %Ŀ�꺯��ֵ����
    for i=1:max(tPhenotype(:,Lind+1)) %��ÿ�ַ������������з����е����ƥ�������
        %��Phenotype������ÿ�е�ǰ2n��Ԫ����n�Խ���ƥ�����������ţ������е�n��ͬ��
        %��������Щƥ������ٵķ���������2n+1:n-1���λ�õ�Ԫ����0,
        %������ǰ��ĸ��ģ���Щ0����1����˶�������ٵķ������������ڲ�������в�õ���ʵ�ǣ�1��1��λ���ϵ�Ԫ�أ���0
        ObjV = ObjV + difference(sub2ind(size(difference), tPhenotype(:,2*i-1), tPhenotype(:,2*i))); %%
    end %�õ�ÿ�ַ��������н���ƥ����߶Ե�ƥ������
    
    ObjV = (T+ObjV)./(tPhenotype(:,Lind+1)); %����ƥ����ۣ�Ŀ�꺯������ʽ�õ�ÿ�ַ�������ƥ�����
    traceM = [traceM,mean(ObjV(find(ObjV~=inf)))]; %��¼�ô���Ⱥ��ƽ������
    traceB = [traceB,min(ObjV)]; %��¼�ô���Ⱥ���������
%--------------------%
    
%--------������Ӧ��--------%
    FitnV = ranking(ObjV); %���������Ӧ��Ϊ2������������������Ӧ�ȣ���i�е�ֵ�ǵ�i�ַ�������Ӧ��
%--------------------%
    
%--------ѡ��--------%
    NewChrIx = rws(FitnV, Nind);
    Chrom = Chrom(NewChrIx, :);
%--------------------%
    
%--------��������--------%
    gen = gen+1;
%--------------------%
end

%%%������ƥ�䷽��%%%
[BFV,BFI] = max(FitnV); %�����Ӧ������(FitnV)�е����ֵ����λ��
BObj = Phenotype(BFI,:);

%%������Ⱥ��������ͼ%%%
figure(figurek),
subplot(2,1,1)
plot(traceM);
title('The trace of the mean object value of population','FontSize',10,'FontAngle','italic');
subplot(2,1,2)
plot(traceB,'.');
title('The trace of the best object value of population','FontSize',10,'FontAngle','italic');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%09.11.05%
