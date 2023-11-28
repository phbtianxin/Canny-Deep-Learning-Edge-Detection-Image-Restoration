function [FittedEdgs,FittedImage] = fitEdgs_f(image,edgPm,BObj,M,N)
%按照最佳匹配方案将各对轮廓线做平滑连接

%其中：
% FittedEdgs --- 各对轮廓线拟合部分的掩膜
% FittedImage --- 将拟合线上色后的图像矩阵
% image --- 一个表示图像的矩阵
% edgPm --- 各轮廓线的特征结构体，在此只使用其中的前两个域变量（与破损区接触的各端点坐标及由其所引的轮廓线坐标序列）
% BObj --- 一个行向量，代表最佳匹配方案
% [M,N] --- 原图像尺寸

%例如：
% [edgPm,w,edgs,M,N,Q] = getEdge_f('planet01.png');
% BObj = optimizeByGA_f(edgPm,M,N);
% image = imread('planet01.png');
% [FittedEdgs,FittedImage] = fitEdgs_f(image,edgPm,BObj,M,N);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FittedEdgs = zeros(M,N);
image = double(image); %转化图像数据类型
for i=1:BObj(end) %对于每对匹配轮廓线
    
    %初始化
    xy1 = [edgPm(BObj(2*i-1)).Curve(:,2),edgPm(BObj(2*i-1)).Curve(:,1)]; %由第端点1所在线段的列序列和行序列分别得到曲线1的x坐标序列和y坐标序列
    xy2 = [edgPm(BObj(2*i)).Curve(:,2),edgPm(BObj(2*i)).Curve(:,1)]; %由第端点2所在线段的列序列和行序列分别得到曲线2的x坐标序列和y坐标序列
    p1 = [edgPm(BObj(2*i-1)).Beginning(2),edgPm(BObj(2*i-1)).Beginning(1)]; %得到一个破损处端点（端点1）的函数坐标
    p2 = [edgPm(BObj(2*i)).Beginning(2),edgPm(BObj(2*i)).Beginning(1)]; %得到另一个破损处端点（端点2）的函数坐标
    
    %%绘制离散点曲线
    figure(i)
    subplot(2,2,1)
    plot(xy1(:,1),xy1(:,2),'.',xy2(:,1),xy2(:,2),'.')
    title('Original discrete curve','FontSize',10,'FontAngle','italic');
    axis equal %使图像的横纵坐标的显示步长一致
    
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
    position1 = find(fix(10*dM5XY1)==0); %查找上述差分值过小的位置
    if length(position1)~=0
        newXY1 = newXY1(1:position1(1)+4,:); %只保留从端点到第一个横坐标5点均值变化过小的点
    end
    
    %砍掉曲线2的重复映射部分
    mean5XY2 = (newXY2(5:end,1)+newXY2(4:end-1,1)+newXY2(3:end-2,1)+newXY2(2:end-3,1)+newXY2(1:end-4,1))/5; %得到第一段曲线x坐标的5点均线
    dM5XY2 = mean5XY2(1:end-1) - mean5XY2(2:end); %做差分。得到每点横坐标均值与下一点横坐标均值的差值
    %position2 = find(dM5XY2(1)*dM5XY2<0); %查找上述差分值符号改变的位置
    position2 = find(fix(10*dM5XY2)==0); %查找上述差分值过小的位置
    if length(position2)~=0
        newXY2 = newXY2(1:position2(1)+4,:); %只保留从端点到第一个横坐标5点均值变化过小的点
    end
     
    %%绘制角度归一化及去重影后的离散点曲线
    figure(i)
    subplot(2,2,2)
    plot(newXY1(:,1),newXY1(:,2),'.',newXY2(:,1),newXY2(:,2),'.')
    title('Unitary discrete curve','FontSize',10,'FontAngle','italic');
    axis equal %使图像的横纵坐标的显示步长一致
    
    %拟合
    f = fit([newXY1(:,1);newXY2(:,1)],[newXY1(:,2);newXY2(:,2)],'poly5'); 
    X = newp1(1):(newp2(1)-newp1(1))/(10*abs(newp2(1)-newp1(1))):newp2(1); %将x坐标也取成小数，从而使拟合出的曲线是一个整个连接分量
    Y = f(X); %求出对应的y坐标

    %%绘制拟合部分的曲线
    figure(i)
    subplot(2,2,3)
    plot(newXY1(:,1),newXY1(:,2),'.',newXY2(:,1),newXY2(:,2),'.')
    hold on
    plot(X,Y,'r.')
    title('Fitted curve','FontSize',10,'FontAngle','italic');
    axis equal %使图像的横纵坐标的显示步长一致
     
    %角度还原
    oXY = [X',Y]*[cos(angle),sin(angle);-sin(angle),cos(angle)]; %得到角度还原后的拟合曲线坐标序列
    oXY = round(oXY); %取整
     
    %%绘制角度还原后的拟合曲线
    figure(i)
    subplot(2,2,4)
    plot(xy1(:,1),xy1(:,2),'.',xy2(:,1),xy2(:,2),'.')
    hold on
    plot(oXY(:,1),oXY(:,2),'r.')
    title('Recovered fitted curve','FontSize',10,'FontAngle','italic');
    axis equal %使图像的横纵坐标的显示步长一致
    
     
    %构造此拟合线在图像中的掩膜
    FittedEdgs(sub2ind(size(FittedEdgs),oXY(:,2),oXY(:,1))) = 1; %将拟合曲线所经位置标1
    
    %在原图上为拟合线上色
    reducedXY = [oXY(find((sum((abs(oXY(1:end-1,:)-oXY(2:end,:)))'))'~=0),:);oXY(end,:)]; %%
    %去掉原拟合曲线坐标序列中的重复坐标。...
    %先差分，找到所有与下一坐标不同的坐标序号；原序列中，这些序号位置上的坐标都是互相不重复的；再添上原序列中的最后一个坐标。
    for l = 1:size(image,3)
        tstep = (image(reducedXY(1,2),reducedXY(1,1),l) - image(reducedXY(end,2),reducedXY(end,1),l))...
            /(size(reducedXY,1)-1); %得到此次色值变化的临时步长
        for j = 2:size(reducedXY,1)-1
            image(reducedXY(j,2),reducedXY(j,1),l) = image(reducedXY(j-1,2),reducedXY(j-1,1),l) - tstep; %延拟合线逐点上色
        end
    end
    
      
    %在拟合线与原线连接处进行预防性插值
    %FittedEdgs(p1(2):(Y(1)-p1(2)+1)/abs(Y(1)-p1(2)+1):Y(1),X(1)) = 1;
    %FittedEdgs(p2(2):(Y(end)-p2(2)+1)/abs(Y(end)-p2(2)+1):Y(end),X(end)) = 1;
    %[此处两端坐标相减后再加1是为了防止坐标重合时出现除以0的情况]
end
FittedEdgs = FittedEdgs==1; %2值化
FittedImage = uint8(image); %转换回图像矩阵的数据类型

%%绘制最终的拟合线掩膜图
%imview(FittedEdgs)

%%绘制含拟合线的照片图像 imview(uint8(image))
imshow(FittedImage)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     