function newM = parallelMove(oriM,mv)

% a = ones(5);
% a = cat(3,a,2*ones(5));
% v = [2 -3];
% b = parallelMove(a,v);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[M,N,L] = size(oriM);

if mv(1)>0 
    oriM(end:-1:end-mv(1)+1,:,:) = []; 
    oriM = [zeros(mv(1),N,L);oriM]; 
end

if mv(1)<0 
    oriM(1:abs(mv(1)),:,:) = []; 
    oriM = [oriM;zeros(abs(mv(1)),N,L)]; 
end

if mv(2)>0 
    oriM(:,end:-1:end-mv(2)+1,:) = []; 
    oriM = [zeros(M,mv(2),L),oriM]; 
end

if mv(2)<0 
    oriM(:,1:abs(mv(2)),:) = []; 
    oriM = [oriM,zeros(M,abs(mv(2)),L)]; 
end

newM = oriM;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

