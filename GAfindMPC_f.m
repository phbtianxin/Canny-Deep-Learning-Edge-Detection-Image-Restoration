function [mpc,Upatch,Dpatch,tp] = GAfindMPC_f(w,edgs,FittedEdgs)
%查找出当前具有最大填充优先级的填充块中心坐标序列
%既要保证找出的这些中心具有最大优先级；
%又要保证任意两个中心之间的横向或是纵向距离相差大于2*tp；
%还要保证每次找出的中心尽可能的多

%前提假设
%1、破损区距图像的边界足够远

%其中：
% mpc --- 一个用来表示当前具有最大填充优先级的填充快中心坐标序列的m行2列矩阵
% Upatch --- 一个m层矩阵，每一层是一个最优填充块中可用部分的掩膜
% Dpatch --- 一个m层矩阵，每一层是一个最优填充块中需填充部分的掩膜
% tp --- 填充块中心到一条边的距离
% w --- 破损区掩膜
% edgs --- 原相关轮廓线掩膜
% FittedEdgs --- 拟合上的轮廓线掩膜

%示例：
% [edgPm,w,edgs,M,N,Q] = getEdge_f('planet01.png');
% BObj = optimizeByGA_f(edgPm,M,N);
% FittedEdgs = fitEdgs_f(edgPm,BObj,M,N);
% [mpc,Upatch,Dpatch,tp] = GAfindMPC_f(w,edgs,FittedEdgs);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[M,N] = size(w);

%tp = 8; %填充块中心到一条边的距离
tp=4;

candidate = bwperim(w,8); %取破损区的最外层像素作为候选中心
candidate = candidate&(~FittedEdgs); %拟合轮廓线上的像素不能被作为候选中心
[canSub(:,1),canSub(:,2)]=find(candidate==1); %获得各候选中心的坐标

%%绘制候选像素图
%imview(candidate)

mpc = []; %初始化具有最大填充优先级的填充块中心坐标序列
mPRI = -1; %初始化最大填充优先级的值。设为-1是为确保Upatch一定能被初始化。
mpcW = zeros(M,N); %初始化一个选中填充块中心的分布图（一个中心一个亮点）
for i=1:size(canSub,1)
    patch = ones(M,N); %初始化一个填充块掩膜
    patch(canSub(i,1)-tp:canSub(i,1)+tp,canSub(i,2)-tp:canSub(i,2)+tp) = 0; %将第i个候选填充块上的像素置0
    patch = patch==1; %2值化
    patch = patch|w|edgs|FittedEdgs; %掩出该填充块中除轮廓线外的各完好区域
    patch(canSub(i,1),canSub(i,2)) = 0; %添上填充块中心
    patch = ~patch; %反白，使各完好区域的值为1
    patch = bwlabel(patch,4); %将各完好区域划分为各连接分量
    patch = patch==patch(canSub(i,1),canSub(i,2)); %将与填充块中心所属的连接分量置1
    patch(canSub(i,1),canSub(i,2)) = 0; %摸去填充块中心
    PRI = sum(sum(patch)); %该填充块的优先级由上述掩膜中的1值像素的个数决定
    
    if PRI>mPRI %如果此填充块的优先级比当前最大填充优先级大
        mPRI = PRI; %则更新当前最大填充优先级
        mpc = canSub(i,:); %重置mpc
        mpcW = zeros(M,N); %重置选中填充块中心的分布图
        mpcW(sub2ind(size(mpcW),mpc(:,1),mpc(:,2))) = 1; %重置选中填充块中心的分布图
        Upatch = []; %重置Upatch
        Upatch(:,:,1) = patch; %重置Upatch
    elseif (PRI==mPRI)&(sum(sum(mpcW(canSub(i,1)-2*tp:canSub(i,1)+2*tp,canSub(i,2)-2*tp:canSub(i,2)+2*tp)))==0) %%
        %当此填充块具有当前最大优先级且不与其他现存填充块有重叠部分
        mpc = [mpc;canSub(i,:)]; %将该填充块中心坐标添加到最优填充块中心坐标序列中
        mpcW(canSub(i,1),canSub(i,2)) = 1; %在现存选中填充块中心分布图上添加该中心点
        Upatch = cat(3,Upatch,patch); %将该填充块中可用部分的掩膜矩阵添加到Upatch中
    end
end
Upatch = Upatch==1; %2值化

%获得各填充块中需填充部分的掩膜
for i=1:size(mpc,1)
    tDpatch = zeros(M,N); %初始化
    tDpatch(mpc(i,1)-tp:mpc(i,1)+tp,mpc(i,2)-tp:mpc(i,2)+tp) = 1; %构造该填充块的整体掩膜
    tDpatch = tDpatch==1; %2值化
    tDpatch = tDpatch&(~FittedEdgs); %用拟合出的轮廓线将该填充块分割成各连接分量（如果有线经过此填充块）
    tDpatch = tDpatch&w; %掩出该填充块中的需填充部分
    tDpatch = bwlabel(tDpatch,4); %划分为各连接分量
    tDpatch = tDpatch==tDpatch(mpc(i,1),mpc(i,2)); %选择该填充中心所在的连接分量
    Dpatch(:,:,i) = tDpatch; %记录
end

%%绘制各填充块中可用部分
%for i=1:size(Upatch,3)
%    imview(Upatch(:,:,i))
%    imview(Dpatch(:,:,i))
%end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    