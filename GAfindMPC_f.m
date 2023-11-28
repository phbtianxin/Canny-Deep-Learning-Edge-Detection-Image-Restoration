function [mpc,Upatch,Dpatch,tp] = GAfindMPC_f(w,edgs,FittedEdgs)

% [edgPm,w,edgs,M,N,Q] = getEdge_f('planet01.png');
% BObj = optimizeByGA_f(edgPm,M,N);
% FittedEdgs = fitEdgs_f(edgPm,BObj,M,N);
% [mpc,Upatch,Dpatch,tp] = GAfindMPC_f(w,edgs,FittedEdgs);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[M,N] = size(w);

tp=4;

candidate = bwperim(w,8); 
candidate = candidate&(~FittedEdgs); 
[canSub(:,1),canSub(:,2)]=find(candidate==1); 

mpc = []; 
mPRI = -1; 
mpcW = zeros(M,N); 
for i=1:size(canSub,1)
    patch = ones(M,N);
    patch(canSub(i,1)-tp:canSub(i,1)+tp,canSub(i,2)-tp:canSub(i,2)+tp) = 0; 
    patch = patch==1; 
    patch = patch|w|edgs|FittedEdgs; 
    patch(canSub(i,1),canSub(i,2)) = 0; 
    patch = ~patch; 
    patch = bwlabel(patch,4); 
    patch = patch==patch(canSub(i,1),canSub(i,2)); 
    patch(canSub(i,1),canSub(i,2)) = 0; 
    PRI = sum(sum(patch)); 
    
    if PRI>mPRI 
        mPRI = PRI; 
        mpc = canSub(i,:); 
        mpcW = zeros(M,N); 
        mpcW(sub2ind(size(mpcW),mpc(:,1),mpc(:,2))) = 1; 
        Upatch = []; 
        Upatch(:,:,1) = patch; 
    elseif (PRI==mPRI)&(sum(sum(mpcW(canSub(i,1)-2*tp:canSub(i,1)+2*tp,canSub(i,2)-2*tp:canSub(i,2)+2*tp)))==0) 
        mpc = [mpc;canSub(i,:)]; 
        mpcW(canSub(i,1),canSub(i,2)) = 1;
        Upatch = cat(3,Upatch,patch); 
    end
end
Upatch = Upatch==1; 

for i=1:size(mpc,1)
    tDpatch = zeros(M,N);
    tDpatch(mpc(i,1)-tp:mpc(i,1)+tp,mpc(i,2)-tp:mpc(i,2)+tp) = 1; 
    tDpatch = tDpatch==1; 
    tDpatch = tDpatch&(~FittedEdgs); 
    tDpatch = tDpatch&w; 
    tDpatch = bwlabel(tDpatch,4);
    tDpatch = tDpatch==tDpatch(mpc(i,1),mpc(i,2)); 
    Dpatch(:,:,i) = tDpatch; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    