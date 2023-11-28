function BObj = optimizeByGA_f(edgPm,M,N)
%根据参与配对的各轮廓线的参数通过GA算法计算最佳匹配方案

%前提假设：
% 1、一定有至少一对相关轮廓线
% 2、认为只有从破损区外向破损区内成汇聚趋势的两条曲线才有匹配资格

%其中，
% BObj -- 一个m+1维行向量（m是相关轮廓线条数）。从第一个元素开始，每相邻两元素是两个要进行配对的轮廓线序号，最后一个元素是配对的对数；空位用0填充。
% edgPm -- 一个Lind×N矩阵，表示Lind条参加配对轮廓线每条线的N个归一化匹配参数。
% M,N -- 图像的尺寸。

%例如：
% [edgPm,w,edgs,M,N,Q] = getEdge_f('planet01.png');
% BObj = optimizeByGA_f(edgPm,M,N);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%参数初始化%%%
Nind = 20; %设置种群个体数
MaxGen = 100; %设置最大遗传代数
T = 1.2; %设置匹配门限
Lind = size(edgPm,2); %获得组成个体的参数个数，即参加匹配的轮廓线条数
%BaseV = [Lind:-1:2,floor(Lind/2)]; %创建基因型向量
BaseV = [Lind:-1:1,floor(Lind/2)]; %创建基因型向量
Chrom = crtbp(Nind,BaseV); %获得初始种群
Pm = 0.7/Lind; %设置变异概率
XOVR = 0.3; %设置交叉概率

traceM = []; %初始化一个种群平均性能跟踪向量
traceB = []; %初始化一个种群最佳性能跟踪向量

%%%计算出每两条轮廓线的匹配差异%%%
difference = zeros(Lind); %初始化一个Lind（轮廓线的条数）阶方阵用来记录每两条轮廓线之间的匹配差异
    
%计算拟合曲线曲折程度
figurek = 1; %初始化绘图序号
Tdifference = 0*difference; %初始化一个临时记录矩阵
for i=1:Lind-1
    for j=i+1:Lind %只需用到difference阵的右上半角即可,不包括主对角线
        
        %初始化
        xy1 = [edgPm(i).Curve(:,2),edgPm(i).Curve(:,1)]; %由第端点1所在线段的列序列和行序列分别得到曲线1的x坐标序列和y坐标序列
        xy2 = [edgPm(j).Curve(:,2),edgPm(j).Curve(:,1)]; %由第端点2所在线段的列序列和行序列分别得到曲线2的x坐标序列和y坐标序列
        p1 = [edgPm(i).Beginning(2),edgPm(i).Beginning(1)]; %得到一个破损处端点（端点1）的函数坐标
        p2 = [edgPm(j).Beginning(2),edgPm(j).Beginning(1)]; %得到另一个破损处端点（端点2）的函数坐标
        
        %%绘制离散点曲线
        %figure(figurek)
        %subplot(2,2,1)
        %plot(xy1(:,1),xy1(:,2),'.',xy2(:,1),xy2(:,2),'.')
        %title('Original discrete curve','FontSize',10,'FontAngle','italic');
        %axis equal %使图像的横纵坐标的显示步长一致
        
        %角度归一化
        v = [p2(1)-p1(1),p2(2)-p1(2)]; %两端点间向量
        angle = atan(v(2)/v(1)); %两端点间向量的相角
        newXY1 = xy1*[cos(angle),-sin(angle);sin(angle),cos(angle)]; %按上述相角旋转曲线1的坐标序列
        newXY2 = xy2*[cos(angle),-sin(angle);sin(angle),cos(angle)]; %按上述相角旋转曲线2的坐标序列
        newp1 = newXY1(1,:); %按相角旋转的端点1坐标
        newp2 = newXY2(1,:); %按相角旋转的端点2坐标
        
        %砍掉曲线1的重复映射部分
        mean5XY1 = (newXY1(5:end,1)+newXY1(4:end-1,1)+newXY1(3:end-2,1)+newXY1(2:end-3,1)+newXY1(1:end-4,1))/5; %得到第一段曲线x坐标的5点均线
        dM5XY1 = mean5XY1(1:end-1) - mean5XY1(2:end); %做差分。得到每点横坐标均值与下一点横坐标均值的差值
        %position1 = find(dM5XY1(1)*dM5XY1<0); %查找上述差分值符号改变的位置
        position1 = find(round(dM5XY1)==0); %查找上述差分值过小的位置
        if length(position1)~=0
            newXY1 = newXY1(1:position1(1)+4,:); %只保留从端点到第一个横坐标5点均值变化过小的点
        end
        
        %砍掉曲线2的重复映射部分
        mean5XY2 = (newXY2(5:end,1)+newXY2(4:end-1,1)+newXY2(3:end-2,1)+newXY2(2:end-3,1)+newXY2(1:end-4,1))/5; %得到第一段曲线x坐标的5点均线
        dM5XY2 = mean5XY2(1:end-1) - mean5XY2(2:end); %做差分。得到每点横坐标均值与下一点横坐标均值的差值
        %position2 = find(dM5XY2(1)*dM5XY2<0); %查找上述差分值符号改变的位置
        position2 = find(round(dM5XY2)==0); %查找上述差分值过小的位置
        if length(position2)~=0
            newXY2 = newXY2(1:position2(1)+4,:); %只保留从端点到第一个横坐标5点均值变化过小的点
        end
        
        %%绘制角度归一化及去重影后的离散点曲线
        %figure(figurek)
        %subplot(2,2,2)
        %plot(newXY1(:,1),newXY1(:,2),'.',newXY2(:,1),newXY2(:,2),'.')
        %title('Unitary discrete curve','FontSize',10,'FontAngle','italic');
        %axis equal %使图像的横纵坐标的显示步长一致
        
        %%%按走向排序线上点坐标后启用
        if ~(((newXY1(end,1)<newXY1(1,1))&&(newXY1(1,1)<newXY2(1,1))&&(newXY2(1,1)<newXY2(end,1)))||...
                ((newXY1(end,1)>newXY1(1,1))&&(newXY1(1,1)>newXY2(1,1))&&(newXY2(1,1)>newXY2(end,1)))) %如果两条线从破损区向外不发散
            Tdifference(i,j) = inf; %则此两条线不适合匹配，从而匹配差异为inf
            continue
        end
        
        %拟合
        %%%考虑自动选择拟合曲线的次数
        f = fit([newXY1(:,1);newXY2(:,1)],[newXY1(:,2);newXY2(:,2)],'poly5');

        %根据其拟合部分曲率的平均变化率求曲线光滑度
        clear cc;
        cc(1) = f.p1; %第1项系数
        cc(2) = f.p2; %第2项系数
        cc(3) = f.p3; %第3项系数
        cc(4) = f.p4; %第4项系数
        cc(5) = f.p5; %第5项系数
        cc(6) = f.p6; %第6项系数
%         cc(7) = f.p7; %第7项系数
%         cc(8) = f.p8; %第8项系数
%         cc(9) = f.p9; %第9项系数
%         cc(10) = f.p10; %第10项系数
        x=sym('x'); %将x定义为一个符号
        cf = cc*[x^5;x^4;x^3;x^2;x;1]; %拟合曲线的符号函数
        cK = abs(diff(diff(cf)))/(1+(diff(cf))^2)^(3/2); %曲线的曲率符号函数
        dcK = diff(cK); %曲线曲率的导函数（符号）
        dx = round(newp1(1)):(newp2(1)-newp1(1))/abs(newp2(1)-newp1(1)):round(newp2(1)); %采样点
        SdcK = abs(subs(dcK,x,dx)); %各采样点处的曲率变化率
        difference(i,j) = sum(SdcK)/length(dx); %整个拟合部分的平均曲率变化率
        
        %%绘制拟合部分的曲线
        %figure(figurek)
        %subplot(2,2,3)
        %plot(newXY1(:,1),newXY1(:,2),'.',newXY2(:,1),newXY2(:,2),'.')
        %hold on
        %plot(dx,(f(dx))','r.')
        %title('Fitted curve','FontSize',10,'FontAngle','italic');
        %axis equal %使图像的横纵坐标的显示步长一致
        
        %%绘制拟合部分的曲率变化率函数图像
        %ScK = subs(cK,x,dx);
        %figure(figurek)
        %subplot(2,2,4)
        %plot(dx,SdcK)
        %title('The derivative of curvature function','FontSize',10,'FontAngle','italic');
        %axis equal %使图像的横纵坐标的显示步长一致
        
        figurek = figurek + 1;
    end
end
difference = difference./(sum(sum(difference))/(sum(sum(difference>0)))); %归一化（即除以该类差异的平均值）
difference = difference + Tdifference;
    
%计算平均色值差异程度
Tdifference = 0*difference; %初始化一个临时记录矩阵
for i=1:Lind
    for j=i+1:Lind %只需用到difference阵的右上半角即可
        Tdifference(i,j) = abs(edgPm(i).MeanColor(1)-edgPm(j).MeanColor(2))+abs(edgPm(i).MeanColor(2)-edgPm(j).MeanColor(1));
    end
end
Tdifference = Tdifference./(sum(sum(Tdifference))/(sum(sum(Tdifference>0))));%归一化（即除以该类差异的平均值）
difference = difference + Tdifference; %求该匹配的各种差异总和
difference = difference + difference'; %将另半角矩阵做个镜像

%%%进行遗传优化%%%
gen = 0; %设置初始代数
while gen < MaxGen
    
%--------变异--------%
    Chrom = mut(Chrom,Pm,BaseV);
%--------------------%
    
%--------交叉--------%
    Chrom = xovsp(Chrom, XOVR);
%--------------------%
    
%--------解码成表现型--------%
    Chrom = Chrom+1; %将基因型以1作为起点
    Phenotype = 0*Chrom; %构造表现型矩阵，除最后一列是配对数，其余每个非0元素表示参加配对的轮廓线序号。从第一个元素，每相邻两元素配对。
    for i=1:Nind
        tem = 1:Lind; %一个用于放置线序号的列表
        for j=1:2*Chrom(i,end) %只需获得要参加配对的线编号即可
            Phenotype(i,j) = tem(Chrom(i,j)); %此处轮廓线的序号是现存序号列表中的Chrom(i,j)位置上的值
            tem(Chrom(i,j)) = []; %已被选走的线序号要从线序号列表中移除
        end
    end
    %Phenotype = [Phenotype,Chrom(:,Lind)]; %加上每种方案需配对的轮廓线对数
    %avoid multi column zero
    Phenotype(:,end)=Chrom(:,Lind+1);
    Chrom = Chrom-1; %将基因型还原为以0作为起点
%--------------------%
    
%--------计算每种方案的匹配代价（目标函数值）--------%
    tPhenotype = Phenotype + (Phenotype==0); %%
    %将表现型矩阵所有的0元素都变成1，从而把查找0号与0号线（即不配对的空位）的匹配差异变成查找1号与1号的匹配差异（会得到0）
    ObjV = zeros(Nind,1); %目标函数值向量
    for i=1:max(tPhenotype(:,Lind+1)) %对每种方案都查找所有方案中的最大匹配对数，
        %在Phenotype矩阵中每行的前2n个元素是n对进行匹配的轮廓线序号，但各行的n不同，
        %所以在有些匹配对数少的方案行里在2n+1:n-1这段位置的元素是0,
        %而经过前面的更改，这些0会变成1，因此对于配对少的方案，到后面在差异矩阵中查得的其实是（1，1）位置上的元素，即0
        ObjV = ObjV + difference(sub2ind(size(difference), tPhenotype(:,2*i-1), tPhenotype(:,2*i))); %%
    end %得到每种方案中所有进行匹配的线对的匹配差异和
    
    ObjV = (T+ObjV)./(tPhenotype(:,Lind+1)); %按照匹配代价（目标函数）公式得到每种方案最后的匹配代价
    traceM = [traceM,mean(ObjV(find(ObjV~=inf)))]; %记录该代种群的平均性能
    traceB = [traceB,min(ObjV)]; %记录该代种群的最佳性能
%--------------------%
    
%--------分配适应度--------%
    FitnV = ranking(ObjV); %采用最大适应度为2的线性排序来分配适应度，第i行的值是第i种方案的适应度
%--------------------%
    
%--------选择--------%
    NewChrIx = rws(FitnV, Nind);
    Chrom = Chrom(NewChrIx, :);
%--------------------%
    
%--------代数增加--------%
    gen = gen+1;
%--------------------%
end

%%%获得最佳匹配方案%%%
[BFV,BFI] = max(FitnV); %获得适应度向量(FitnV)中的最大值及其位置
BObj = Phenotype(BFI,:);

%%绘制种群性能趋势图%%%
figure(figurek),
subplot(2,1,1)
plot(traceM);
title('The trace of the mean object value of population','FontSize',10,'FontAngle','italic');
subplot(2,1,2)
plot(traceB,'.');
title('The trace of the best object value of population','FontSize',10,'FontAngle','italic');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%09.11.05%
