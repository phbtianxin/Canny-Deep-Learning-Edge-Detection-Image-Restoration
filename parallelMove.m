function newM = parallelMove(oriM,mv)
%将一个矩阵在原尺寸范围内做平移。移出的部分丢弃，移入的部分用0填充

%其中：
% newM --- 与输入矩阵同尺寸的平移后的矩阵
% oriM --- 要做平移处理的输入矩阵。如果是三维的，则视为由多个同尺寸的2维矩阵罗列在一起同时做平移
% mv --- 一个2元向量，作为平抑向量。即该向量由原位置指向对应的新位置。其第一个元表示行变化，第二个元表示列变化。

%例如:
% a = ones(5);
% a = cat(3,a,2*ones(5));
% v = [2 -3];
% b = parallelMove(a,v);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[M,N,L] = size(oriM);

if mv(1)>0 %此时要下移
    oriM(end:-1:end-mv(1)+1,:,:) = []; %在矩阵底部去掉mv(1)行
    oriM = [zeros(mv(1),N,L);oriM]; %在矩阵顶部添加mv(1)行
end

if mv(1)<0 %此时要上移
    oriM(1:abs(mv(1)),:,:) = []; %在矩阵顶部去掉abs(mv(1))行
    oriM = [oriM;zeros(abs(mv(1)),N,L)]; %在矩阵底部添加abs(mv(1))行
end

if mv(2)>0 %此时要右移
    oriM(:,end:-1:end-mv(2)+1,:) = []; %在矩阵右侧去掉mv(2)列
    oriM = [zeros(M,mv(2),L),oriM]; %在矩阵左侧添加mv(2)列
end

if mv(2)<0 %此时要左移
    oriM(:,1:abs(mv(2)),:) = []; %在矩阵左侧去掉abs(mv(2))列
    oriM = [oriM,zeros(M,abs(mv(2)),L)]; %在矩阵右侧添加abs(mv(2))列
end

newM = oriM;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

