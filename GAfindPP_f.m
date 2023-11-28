function [ppc,Rpatch] = GAfindPP_f(image,w,edgs,FittedEdgs,mpc,Upatch,Dpatch,tp)
%查找出对应于各填充块的最佳匹配块的中心坐标
%使用外部函数：parallelMove()

%前提假设
%1、破损区距图像的边界足够远
%2、在扫描区域内只少存在与填充块尺寸一致的完整匹配块

%其中：
% ppc --- 一个m行2列矩阵，其第i行记录的是对应于第i个填充块的最佳匹配块的中心坐标
% Rpatch --- 一个m层矩阵，每一层是一个最优匹配块中可填充部分的掩膜
% image --- 一个表示图像的矩阵
% w --- 破损区掩膜
% edgs --- 原相关轮廓线掩膜
% FittedEdgs --- 拟合上的轮廓线掩膜
% mpc --- 一个用来表示当前具有最大填充优先级的填充快中心坐标序列的m行2列矩阵
% Upatch --- 一个m层矩阵，每一层是一个最优填充块中可用部分的掩膜
% Dpatch --- 一个m层矩阵，每一层是一个最优填充块中需填充部分的掩膜
% tp --- 填充块中心到填充块任一条边的距离

%示例:
% [edgPm,w,edgs,M,N,Q] = getEdge_f('planet01.png');
% BObj = optimizeByGA_f(edgPm,M,N);
% FittedEdgs = fitEdgs_f(edgPm,BObj,M,N);
% [mpc,Upatch,Dpatch,tp] = GAfindMPC_f(w,edgs,FittedEdgs,5);
% image = imread('planet01.png');
% [ppc,Rpatch] = GAfindPP_f(image,w,edgs,FittedEdgs,mpc,Upatch,Dpatch,tp);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tr = 8; %填充块中心到扫描区任一条边所对应的扫描中心的距离（必须要不tp大）

[M,N,L] = size(image);
uw = w|edgs|FittedEdgs; %构造一个包含破损区及相关轮廓线的掩膜（可用区掩膜）

ppc = []; %初始化最佳匹配块的中心坐标序列
for k = 1 : size(mpc,1) %对于给定的各填充块
    
    %获得对于该填充中心的扫描区域内的可用部分的掩膜
    %前提：相关轮廓线足够分割当前扫描区域
    region = zeros(M,N); %初始化可用扫描区掩膜（即不包括破损区、相关轮廓线、及非同一连接分量）
    region(mpc(k,1)-tr-tp:mpc(k,1)+tr+tp,mpc(k,2)-tr-tp:mpc(k,2)+tr+tp) = 1; %构造该填充块扫描区掩膜
    region = region==1; %2值化
    region = region&(~uw); %在可用区掩膜中掩出该扫描区内的可用部分
    region(mpc(k,1),mpc(k,2)) = 1; %添上该填充块中心
    region = bwlabel(region,4); %将各可用区域划分为各连接分量
    region = region==region(mpc(k,1),mpc(k,2)); %将与该填充块中心所属的连接分量置1
    region(mpc(k,1),mpc(k,2)) = 0; %抹去该填充块中心
    
    UDp = Upatch(:,:,k)|Dpatch(:,:,k); %得到该填充块的可用部分和需填充部分的联合掩膜
    UDp = UDp(mpc(k,1)-tp:mpc(k,1)+tp,mpc(k,2)-tp:mpc(k,2)+tp); %抠出上述掩膜中的填充块部分
    region = imerode(region,UDp); %用该填充块的可用部分和需填充部分一起腐蚀初始扫描区，得到完全适合该填充块的扫描区域
    regionSub = []; %重新初始化
    [regionSub(:,1),regionSub(:,2)] = find(region==1); %提取上述扫描区坐标序列

    minSV = inf; %%初始化最小平均误差平方和
    ppc(k,:) = [0,0]; %初始化对应于第i个填充块的最佳匹配块的中心坐标
    Rpatch(:,:,k) = zeros(M,N); %初始化对应于该填充块的匹配块中可填充部分的掩膜
    for l = 1:size(regionSub,1)
        i = regionSub(l,1); %行坐标
        j = regionSub(l,2); %列坐标
        
        %将该填充块的可用部分和需填充部分同时平移到当前扫描位置
        Upw = parallelMove(Upatch(:,:,k),[i-mpc(k,1),j-mpc(k,2)]); %得到平移该填充块中可用部分到当前扫描位置的掩膜
        Dpw = parallelMove(Dpatch(:,:,k),[i-mpc(k,1),j-mpc(k,2)]); %得到平移该填充块中需填充部分到当前扫描位置的掩膜
        Upw = Upw==1; %2值化
        Dpw = Dpw==1; %2值化
        
        if L==3
            Upw = cat(3,Upw,cat(3,Upw,Upw)); %扩展成与输入图像同维
            oUpw = cat(3,Upatch(:,:,k),cat(3,Upatch(:,:,k),Upatch(:,:,k))); %扩展成与输入图像同维
        else
            oUpw = Upatch(:,:,k); %灰度图像处理
        end
        
        Up = Upw.*double(image); %在输入图像中掩出当前扫描位置中的可用部分
        Up = parallelMove(Up,[mpc(k,1)-i,mpc(k,2)-j]); %将上述掩出的部分再平移回该填充块位置
        oUp = oUpw.*double(image); %在输入图像中掩出当该填充块中的可用部分
        SV = (sum(sum(sum((Up-oUp).^2))))/sum(sum(Upatch(:,:,k))); %计算当前扫描块与填充块中同属可用部分之间的平均误差平方和
        
        if SV<minSV %如果该平均误差平方和小于当前最小平均误差平方和
            minSV = SV; %%更新当前最小平均误差平方和
            Rpatch(:,:,k) = Dpw; %更新对应于该填充块的匹配块中可填充部分的掩膜
            ppc(k,:) = [i,j]; %更新对应于第i个填充块的最佳匹配块的中心坐标
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
    