function [edgPm,w,edgs,M,N,Q] = getEdge_f(image)
%提取含破损区域图像中与破损区域有接触的各条轮廓线序列及相应各描绘子

%其中：
% edgPm --- 一个结构体序列（有几个与破损区接触的轮廓线端点就有几个这样的结构体），每个结构体的域变量包括：
% 1、与破损区相接触的一条轮廓线的端点坐标
% 2、该端点所在条轮廓线的坐标序列
% 3、该部分轮廓线左右两侧的平均亮度
% 4、该部分轮廓线左右两侧的对比度（方差）
% image --- 一个字符串，代表一幅图像
% w --- 图像破损区域的掩膜
% edgs --- 所有相关轮廓线的掩膜
% [M,N,Q] --- 图像的尺寸

%示例：
% [edgPm,w,edgs,M,N,Q] = getEdge_f('planet01.png');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%image = imread(image); %获取图像矩阵
[M,N,Q] = size(image); %获得图像的尺寸
if Q==3
    image = rgb2gray(image); %获得图像的灰度值矩阵
end

w = image==0; %获得图像破损区域的掩膜

%获得描述图像中轮廓的2值矩阵
g = edge(image, 'canny', [0.4, 0.5], 1.4); 
% 'canny'后中括号中的两个值越大，找出的轮廓越强；最后一个值越大，找出的轮廓越平滑
imshow(g);
%将当前破损区的掩膜
w_N = [w(2:M,:);zeros(1,N)]==1; %上挫一行
w_S = [zeros(1,N);w(1:M-1,:)]==1; %下挫一行
w_W = [w(:,2:N),zeros(M,1)]==1; %左挫一列
w_E = [zeros(M,1),w(:,1:N-1)]==1; %右挫一列
w_NW = [w_N(:,2:N),zeros(M,1)]==1; %上挫一行并左挫一列
w_NE = [zeros(M,1),w_N(:,1:N-1)]==1; %上挫一行并右挫一列
w_SW = [w_S(:,2:N),zeros(M,1)]==1; %下挫一行并左挫一列
w_SE = [zeros(M,1),w_S(:,1:N-1)]==1; %下挫一行并右挫一列
bw = w_N|w_S|w_W|w_E|w_NW|w_NE|w_SW|w_SE|w; %将破损区掩膜扩一圈边
newg = (g-bw)==1; %去掉破损区掩膜的轮廓

%将当前破损区的掩膜
w_N = [bw(2:M,:);zeros(1,N)]==1; %上挫一行
w_S = [zeros(1,N);bw(1:M-1,:)]==1; %下挫一行
w_W = [bw(:,2:N),zeros(M,1)]==1; %左挫一列
w_E = [zeros(M,1),bw(:,1:N-1)]==1; %右挫一列
w_NW = [w_N(:,2:N),zeros(M,1)]==1; %上挫一行并左挫一列
w_NE = [zeros(M,1),w_N(:,1:N-1)]==1; %上挫一行并右挫一列
w_SW = [w_S(:,2:N),zeros(M,1)]==1; %下挫一行并左挫一列
w_SE = [zeros(M,1),w_S(:,1:N-1)]==1; %下挫一行并右挫一列
bbw = w_N|w_S|w_W|w_E|w_NW|w_NE|w_SW|w_SE|bw; %将破损区掩膜再扩一圈边

%将当前破损区的掩膜
w_N = [bbw(2:M,:);zeros(1,N)]==1; %上挫一行
w_S = [zeros(1,N);bbw(1:M-1,:)]==1; %下挫一行
w_W = [bbw(:,2:N),zeros(M,1)]==1; %左挫一列
w_E = [zeros(M,1),bbw(:,1:N-1)]==1; %右挫一列
w_NW = [w_N(:,2:N),zeros(M,1)]==1; %上挫一行并左挫一列
w_NE = [zeros(M,1),w_N(:,1:N-1)]==1; %上挫一行并右挫一列
w_SW = [w_S(:,2:N),zeros(M,1)]==1; %下挫一行并左挫一列
w_SE = [zeros(M,1),w_S(:,1:N-1)]==1; %下挫一行并右挫一列
bbbw = w_N|w_S|w_W|w_E|w_NW|w_NE|w_SW|w_SE|bbw; %将破损区掩膜再扩一圈边

%选择与破损区域有接触的轮廓线
tbind = find(bw==1); %找到破损区中第一个像素的坐标
tL = bwlabel(bbw|newg,8); %将与破损区相接触的轮廓线及破损区合并成一个连接分量
tL = tL==tL(tbind(1)); %仅将上述的连接分量的值设成1，其余全是0，从而得到一个只包含破损区及与破损区相接触轮廓线的掩膜
finalg = newg&tL; %用上述掩膜掩出与破损区相接触的轮廓线
tL = bwlabel(finalg,8); %表示各相关轮廓线
for i=1:max(max(tL))
    [tI,tJ] = find(tL==i);
    edgeCurve{i,1} = [tI,tJ];
end %查找每条轮廓线的坐标序列

edgs = zeros(M,N); %初始化相关轮廓线掩膜
for i=1:length(edgeCurve)
    edgs(sub2ind([M,N],edgeCurve{i,1}(:,1),edgeCurve{i,1}(:,2))) = 1;
end %将所有相关轮廓线的像素置1
edgs = edgs==1; %2值化

%%绘制所提取的轮廓线图
%imview(edgs)
%imview(w)
imshow(edgs)

%各条相关轮廓线的各描绘子
d = 5; %描绘子提取框的宽度
k = 1;
AllC = [];
for i=1:length(edgeCurve)
    finalg = zeros(M,N); %初始化第i条轮廓线掩膜
    finalg(sub2ind([M,N], edgeCurve{i,1}(:,1), edgeCurve{i,1}(:,2)))=1; %标出第i条轮廓线
    finalg = finalg==1; %2值化
    
    n = 1;
    for j=1:size(edgeCurve{i,1},1) %对于某轮廓线上的所有像素
        if sum(sum(finalg(edgeCurve{i,1}(j,1)-1:edgeCurve{i,1}(j,1)+1,edgeCurve{i,1}(j,2)-1:edgeCurve{i,1}(j,2)+1)))<3
            Cind(n,:) = edgeCurve{i,1}(j,:); %则该像素即为端点像素
            %若某一像素周围与之相连的像素只有一个
            n = n+1;
        elseif sum(sum(finalg(edgeCurve{i,1}(j,1)-1:edgeCurve{i,1}(j,1)+1,edgeCurve{i,1}(j,2)-1:edgeCurve{i,1}(j,2)+1)))==3
             if (sum(sum(finalg(edgeCurve{i,1}(j,1)-1:edgeCurve{i,1}(j,1),edgeCurve{i,1}(j,2)-1:edgeCurve{i,1}(j,2)+1)))==1)||(sum(sum(finalg(edgeCurve{i,1}(j,1):edgeCurve{i,1}(j,1)+1,edgeCurve{i,1}(j,2)-1:edgeCurve{i,1}(j,2)+1)))==1)
                  %delete middle point of the edgeCurve 
                  if (edgeCurve{i,1}(j,2)-1==1)&&(sum(finalg(edgeCurve{i,1}(j,1),edgeCurve{i,1}(j,2):edgeCurve{i,1}(j,2)+2))<2)
                       Cind(n,:) = edgeCurve{i,1}(j,:); %则该像素即为端点像素
                       n = n+1;
                  elseif (edgeCurve{i,1}(j,2)+1==N)&&(sum(finalg(edgeCurve{i,1}(j,1),edgeCurve{i,1}(j,2)-2:edgeCurve{i,1}(j,2)))<2)
                       Cind(n,:) = edgeCurve{i,1}(j,:); %则该像素即为端点像素
                       n = n+1;
                  elseif (edgeCurve{i,1}(j,2)-1>1)&&(edgeCurve{i,1}(j,2)+1<N)
                      if (sum(finalg(edgeCurve{i,1}(j,1),edgeCurve{i,1}(j,2)-2:edgeCurve{i,1}(j,2)))<2)&&(sum(finalg(edgeCurve{i,1}(j,1),edgeCurve{i,1}(j,2):edgeCurve{i,1}(j,2)+2))<2)
                          Cind(n,:) = edgeCurve{i,1}(j,:); %则该像素即为端点像素
                          n = n+1;
                      else
                          finalg(edgeCurve{i,1}(j,1),edgeCurve{i,1}(j,2))=0;
                      end  
                  end                  
             end
             if (sum(sum(finalg(edgeCurve{i,1}(j,1)-1:edgeCurve{i,1}(j,1)+1,edgeCurve{i,1}(j,2):edgeCurve{i,1}(j,2)+1)))==1)||(sum(sum(finalg(edgeCurve{i,1}(j,1)-1:edgeCurve{i,1}(j,1)+1,edgeCurve{i,1}(j,2)-1:edgeCurve{i,1}(j,2))))==1)
                 %delete middle point of the edgeCurve 
                 if (edgeCurve{i,1}(j,1)-1==1)&&(sum(finalg(edgeCurve{i,1}(j,1):edgeCurve{i,1}(j,1)+2,edgeCurve{i,1}(j,2)))<2)
                     Cind(n,:) = edgeCurve{i,1}(j,:); %则该像素即为端点像素
                     n = n+1;
                 elseif (edgeCurve{i,1}(j,1)+1==M)&&(sum(finalg(edgeCurve{i,1}(j,1)-2:edgeCurve{i,1}(j,1),edgeCurve{i,1}(j,2)))<2)
                     Cind(n,:) = edgeCurve{i,1}(j,:); %则该像素即为端点像素
                     n = n+1;                     
                 elseif (edgeCurve{i,1}(j,1)-1>1)&&(edgeCurve{i,1}(j,1)+1<M)
                     if (sum(finalg(edgeCurve{i,1}(j,1)-2:edgeCurve{i,1}(j,1),edgeCurve{i,1}(j,2)))<2)&&(sum(finalg(edgeCurve{i,1}(j,1):edgeCurve{i,1}(j,1)+2,edgeCurve{i,1}(j,2)))<2)
                         Cind(n,:) = edgeCurve{i,1}(j,:); %则该像素即为端点像素
                         n = n+1;                                         
                     else
                         finalg(edgeCurve{i,1}(j,1),edgeCurve{i,1}(j,2))=0;
                     end
                 end
             end
        end
    end%获得第i条轮廓线两端的坐标

    %%为绘制各相关轮廓线的端点
    %AllC = [AllC;Cind];
    
    if bbbw(Cind(1,1),Cind(1,2)) %如果端点1与破损区域接触
        pw = zeros(M,N); %初始化端点描绘子提取框掩膜
        
        tw = pw; %为获得轮廓线上端点的前一点坐标初始化一个临时掩膜
        tw(Cind(1,1)-1:Cind(1,1)+1,Cind(1,2)-1:Cind(1,2)+1) = 1; %在端点处扩一圈边
        tw(Cind(1,1),Cind(1,2)) = 0; %摸去端点
        tw = tw==1; %2值化
        tw = tw&finalg; %掩出轮廓线上端点的前一点
        [R,C] = find(tw==1); %得到廓线上端点的前一点的坐标
        %Two point method
        if size([R,C],1)>1
            R=R(1);
            C=C(1);
        end  
        
        anticR = Cind(1,1) + round(((R-Cind(1,1))-(C-Cind(1,2)))/2); %得到以端点为中心逆时针转的下一点的行坐标
        anticC = Cind(1,2) + round(((R-Cind(1,1))+(C-Cind(1,2)))/2); %得到以端点为中心逆时针转的下一点的列坐标
        
        pw(Cind(1,1)-d:Cind(1,1)+d,Cind(1,2)-d:Cind(1,2)+d) = 1; %假设任意两条轮廓线之间的距离都大于描绘子提取框的宽度（不大于又如何？）
        pw = pw==1; %得到以轮廓线端点为中心的一个方形窗口掩膜
        
        region = pw&~(bbw|finalg); %标出方形窗口中除轮廓线及一部分破损区外的部分
       
        L = bwlabel(region,4); %得到方形窗口中由轮廓线及一部分破损区所分开的两个连接分量
        pwRight = L==L(anticR,anticC); %获得从轮廓线向破损区看去的右侧部分的提取掩膜
        pwLeft = region&(~pwRight); %获得从轮廓线向破损区看去左侧部分的提取掩膜
        
        %%绘制从i条轮廓线向破损区看去左右两侧部分的提取掩膜
        %imview(pwRight);
        %imview(pwLeft);
        
        edgPm(k).Beginning = [Cind(1,1),Cind(1,2)]; %将第i条轮廓线与破损区的接触端坐标记录到输出结构的Beginning域变量中
        
        %按线形走向对该轮廓线的坐标序列重排序
        finalg(edgPm(k).Beginning(1),edgPm(k).Beginning(2)) = 0; %去掉端点
        clear lineSub

        [lineSub(:,1),lineSub(:,2)] = find(finalg==1); %找到线上所有点坐标
        newCurve = edgPm(k).Beginning; %将端点坐标置于新序列的第一行
        while size(lineSub,1)~=0
            distance = abs(ones(size(lineSub,1),1)*newCurve(end,:)-lineSub); %计算前一点与余下各点的距离
            sumd = distance(:,1).^2 + distance(:,2).^2; %得到横纵坐标差值的平方和
            if find(sumd<=1)
                newCurve = [newCurve;lineSub(find(sumd<=1),:)]; %将找到的下一点坐标添加到心序列中
            end
            if find(sumd>1&sumd<3)
                newCurve = [newCurve;lineSub(find(sumd>1&sumd<3),:)]; %将找到的下一点坐标添加到心序列中            
            end   
            lineSub(find(sumd<3),:) = []; %在旧序列中摸去该点坐标
        end
        
        edgPm(k).Curve = newCurve; %将第i条轮廓线序列记录到输出结构的Curve域变量中
        edgPm(k).MeanColor = [sum(sum(pwLeft.*double(image)))/sum(sum(pwLeft)),sum(sum(pwRight.*double(image)))/sum(sum(pwRight))]; %%
        %将从第i条轮廓线向破损区看去的左右两侧图像的平均亮度记录到输出结构的MeanColor域变量中
        edgPm(k).variance = [sqrt((sum(sum((pwLeft.*(double(image)-edgPm(k).MeanColor(1,1))).^2)))/sum(sum(pwLeft))),...
            sqrt((sum(sum((pwRight.*(double(image)-edgPm(k).MeanColor(1,2))).^2)))/sum(sum(pwRight)))];
        k = k+1;
    end
    
    if bbbw(Cind(2,1),Cind(2,2)) %如果端点2与破损区域接触
        pw = zeros(M,N); %初始化端点描绘子提取框掩膜
        
        tw = pw; %为获得轮廓线上端点的前一点坐标初始化一个临时掩膜
        tw(Cind(2,1)-1:Cind(2,1)+1,Cind(2,2)-1:Cind(2,2)+1) = 1; %在端点处扩一圈边
        tw(Cind(2,1),Cind(2,2)) = 0; %摸去端点
        tw = tw==1; %2值化
        tw = tw&finalg; %掩出轮廓线上端点的前一点
        [R,C] = find(tw==1); %得到廓线上端点的前一点的坐标
        %Two point method
        if size([R,C],1)>1
            R=R(1);
            C=C(1);
        end               
        anticR = Cind(2,1) + round(((R-Cind(2,1))-(C-Cind(2,2)))/2); %得到以端点为中心逆时针转的下一点的行坐标
        anticC = Cind(2,2) + round(((R-Cind(2,1))+(C-Cind(2,2)))/2); %得到以端点为中心逆时针转的下一点的列坐标
        
        pw(Cind(2,1)-d:Cind(2,1)+d,Cind(2,2)-d:Cind(2,2)+d) = 1; %假设任意两条轮廓线之间的距离都大于描绘子提取框的宽度
        pw = pw==1; %得到以轮廓线端点为中心的一个方形窗口掩膜
        
        region = pw&~(bbw|finalg); %标出方形窗口中除轮廓线及一部分破损区外的部分
        
        L = bwlabel(region,4); %得到方形窗口中由轮廓线及一部分破损区所分开的两个连接分量
        pwRight = L==L(anticR,anticC); %获得从轮廓线向破损区看去的右侧部分的提取掩膜
        pwLeft = region&(~pwRight); %获得从轮廓线向破损区看去左侧部分的提取掩膜
        
        %%绘制从i条轮廓线向破损区看去左右两侧部分的提取掩膜
        %imview(pwRight);
        %imview(pwLeft);
        
        edgPm(k).Beginning = [Cind(2,1),Cind(2,2)]; %将第i条轮廓线与破损区的接触端坐标记录到输出结构的Beginning域变量中
        
        %按线形走向对该轮廓线的坐标序列重排序
        finalg(edgPm(k).Beginning(1),edgPm(k).Beginning(2)) = 0; %去掉端点
        clear lineSub
        [lineSub(:,1),lineSub(:,2)] = find(finalg==1); %找到线上所有点坐标
        newCurve = edgPm(k).Beginning; %将端点坐标置于新序列的第一行
        while size(lineSub,1)~=0
            distance = abs(ones(size(lineSub,1),1)*newCurve(end,:)-lineSub); %计算前一点与余下各点的距离
            sumd = distance(:,1).^2 + distance(:,2).^2; %得到横纵坐标差值的平方和
            %write oder rom small to large
            if find(sumd<=1)
                newCurve = [newCurve;lineSub(find(sumd<=1),:)]; %将找到的下一点坐标添加到心序列中
            end
            if find(sumd>1&sumd<3)
                newCurve = [newCurve;lineSub(find(sumd>1&sumd<3),:)]; %将找到的下一点坐标添加到心序列中            
            end            
            lineSub(find(sumd<3),:) = []; %在旧序列中摸去该点坐标
        end
        
        edgPm(k).Curve = newCurve; %将第i条轮廓线序列记录到输出结构的Curve域变量中
        edgPm(k).MeanColor = [sum(sum(pwLeft.*double(image)))/sum(sum(pwLeft)),sum(sum(pwRight.*double(image)))/sum(sum(pwRight))]; %%
        %将从第i条轮廓线向破损区看去的左右两侧图像的平均亮度记录到输出结构的MeanColor域变量中
        edgPm(k).variance = [sqrt((sum(sum((pwLeft.*(double(image)-edgPm(k).MeanColor(1,1))).^2)))/sum(sum(pwLeft))),...
            sqrt((sum(sum((pwRight.*(double(image)-edgPm(k).MeanColor(1,2))).^2)))/sum(sum(pwRight)))];
        k = k+1;
    end
end

%%绘制个相关轮廓线的端点
%drawp = zeros(M,N); %初始化
%drawp(sub2ind([M,N], AllC(:,1), AllC(:,2))) = 1; %将端点标1
%drawp = drawp==1; %2值化
%imview(drawp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


