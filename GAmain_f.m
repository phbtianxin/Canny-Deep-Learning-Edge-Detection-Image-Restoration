function RImage = GAmain_f(Image)
%对一幅图画的破损区（由纯黑标示）进行基于遗传算法的修复

%其中：
% RImage --- 一个表示图像的矩阵
% Image --- 一个字符串，代表一幅图像

%示例：
% RImage = GAmain_f('planet01.png');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%image = imread(Image);

%%绘制原图像
%imview(image)
image=Image;
imshow(image);%function imview has been removed 

%提取含破损区域图像中与破损区域有接触的各条轮廓线序列及相应各描绘子
[edgPm,w,edgs,M,N,Q] = getEdge_f(Image);

%根据参与配对的各轮廓线的参数通过GA算法计算最佳匹配方案
BObj = optimizeByGA_f(edgPm,M,N);

%按照最佳匹配方案将各对轮廓线做平滑连接
[FittedEdgs,FittedImage] = fitEdgs_f(image,edgPm,BObj,M,N);

RImage = double(FittedImage); %初始化输出图像矩阵
%avoid endless loop
Record=[(sum(sum(w&~FittedEdgs))),1];
while sum(sum(w&~FittedEdgs))~=0 %只要除拟合的轮廓线外还存在破损区
    
    %查找出当前具有最大填充优先级的填充块中心坐标序列
    [mpc,Upatch,Dpatch,tp] = GAfindMPC_f(w,edgs,FittedEdgs);
    
    %查找出对应于各填充块的最佳匹配块的中心坐标
    [ppc,Rpatch] = GAfindPP_f(RImage,w,edgs,FittedEdgs,mpc,Upatch,Dpatch,tp);
    
    for i=1:size(ppc,1)
        Rpw = Rpatch(:,:,i); %得到目标匹配块中的可填充部分的掩膜
        
        if Q==3 %如果是彩色图像
            Rpw = cat(3,Rpw,cat(3,Rpw,Rpw)); %扩展成与输入图像同维
        end
        
        Rp = Rpw.*RImage; %在输入图像中掩出用于填充该破损区的部分
        Rp = parallelMove(Rp,[mpc(i,1)-ppc(i,1),mpc(i,2)-ppc(i,2)]); %将上述掩出的部分再平移回需填充的位置
        RImage = RImage + Rp; %填充图像
        w = w&(~(Rp(:,:,1)>0));%填充破损区掩膜
    end
    %avoid endless loop
    if Record(end,2)>10
        break;
    else
        if sum(sum(w&~FittedEdgs))==Record(end,1)
            Record(end,2)=Record(end,2)+1;
        else
            Record=[Record;[(sum(sum(w&~FittedEdgs))),1]];
        end
    end
end
RImage = uint8(RImage); %还原成图像的数据格式

%%绘制修复后的图像 imview(uint8(RImage))
%imview(RImage)
imshow(RImage);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        