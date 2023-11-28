function ObjV = objFun_f(Chrom,edgPm)
%根据参与配对的各轮廓线的参数计算几种配对方案的匹配代价。

%其中，
% ObjV -- 一个列向量，其第i行值对应第i种配对方案的匹配代价；
% Chrom -- 一个Nind×Lind矩阵，Nind表示染色体个数，即配对方案个数；
%          Lind表示染色体中的基因个数长度，此处每个基因用一个整数表
%          示，其中前Lind-1个代表这Lind条轮廓线的一种排列顺序，最后
%          一个参数代表在这种排序中，前多少对线是被两两配对的；
% edgPm -- 一个Lind×N矩阵，表示Lind条参加配对轮廓线每条线的N个归一化匹配参数。

%例如：
% N = 5;
% Chrom = crtbp(7, [N N-1 N-2 N-3 floor(N/2)]);
% edgPm = rand(5,4);
% ObjV = objFun_f(Chrom,edgPm);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T = 0.3; %设置匹配门限
[Nind,Lind] = size(Chrom);

%解码成表现型
Chrom = Chrom+1;
Phenotype = 0*Chrom; 
for i=1:Nind
    tem = 1:Lind; %一个用于放置线序号的列表
    for j=1:2*Chrom(i,Lind) %只需获得要参加配对的线编号即可
        Phenotype(i,j) = tem(Chrom(i,j)); %此处轮廓线的序号是现存序号列表中的Chrom(i,j)位置上的值
        tem(Chrom(i,j)) = []; %已被选走的线序号要从线序号列表中移除
    end
end
Nc = Chrom(:,Lind); %需配对数向量
        

%构造加权矩阵来对每组进行配对的两条线计算其参数的匹配差异
ObjV = zeros(Nind,1);
k = 0;
while(sum(Nc)~=0)
    R = 0*Phenotype; %除去进行配对的两条线的参数外，其他线的参数的权值都为0
    for i=1:Nind
        if Nc(i,1) ~= 0
            R(i,Phenotype(i,k+1))=1; %配对的第一条线的所有参数被加以权值1
            R(i,Phenotype(i,k+2))=-1; %配对的第二条线的所有参数被加以权值-1
            Nc(i,1) = Nc(i,1)-1;
        end
    end
    ObjV = ObjV+(sum((abs(R*edgPm))'))'; %对于每种匹配方案，总匹配差异等于前一对匹配差异加新一对的匹配差异
    k = k+2;
end

%对于每种匹配方案，用差异和加匹配门限再除匹配对数，即得目标函数
ObjV = (T+ObjV)./(Chrom(:,Lind));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%09.11.04%

