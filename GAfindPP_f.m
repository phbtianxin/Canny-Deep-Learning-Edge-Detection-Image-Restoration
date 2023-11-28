function [ppc,Rpatch] = GAfindPP_f(image,w,edgs,FittedEdgs,mpc,Upatch,Dpatch,tp)
% [edgPm,w,edgs,M,N,Q] = getEdge_f('planet01.png');
% BObj = optimizeByGA_f(edgPm,M,N);
% FittedEdgs = fitEdgs_f(edgPm,BObj,M,N);
% [mpc,Upatch,Dpatch,tp] = GAfindMPC_f(w,edgs,FittedEdgs,5);
% image = imread('planet01.png');
% [ppc,Rpatch] = GAfindPP_f(image,w,edgs,FittedEdgs,mpc,Upatch,Dpatch,tp);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tr = 8; 

[M,N,L] = size(image);
uw = w|edgs|FittedEdgs; 

ppc = [];
for k = 1 : size(mpc,1) 
    
    region = zeros(M,N); 
    region(mpc(k,1)-tr-tp:mpc(k,1)+tr+tp,mpc(k,2)-tr-tp:mpc(k,2)+tr+tp) = 1; 
    region = region==1; 
    region = region&(~uw); 
    region(mpc(k,1),mpc(k,2)) = 1; 
    region = bwlabel(region,4); 
    region = region==region(mpc(k,1),mpc(k,2)); 
    region(mpc(k,1),mpc(k,2)) = 0; 
    
    UDp = Upatch(:,:,k)|Dpatch(:,:,k); 
    UDp = UDp(mpc(k,1)-tp:mpc(k,1)+tp,mpc(k,2)-tp:mpc(k,2)+tp); 
    region = imerode(region,UDp); 
    regionSub = []; 
    [regionSub(:,1),regionSub(:,2)] = find(region==1); 

    minSV = inf;
    ppc(k,:) = [0,0]; 
    Rpatch(:,:,k) = zeros(M,N); 
    for l = 1:size(regionSub,1)
        i = regionSub(l,1); 
        j = regionSub(l,2);
        
        Upw = parallelMove(Upatch(:,:,k),[i-mpc(k,1),j-mpc(k,2)]); 
        Dpw = parallelMove(Dpatch(:,:,k),[i-mpc(k,1),j-mpc(k,2)]); 
        Upw = Upw==1; 
        Dpw = Dpw==1; 
        
        if L==3
            Upw = cat(3,Upw,cat(3,Upw,Upw)); 
            oUpw = cat(3,Upatch(:,:,k),cat(3,Upatch(:,:,k),Upatch(:,:,k))); 
        else
            oUpw = Upatch(:,:,k); 
        end
        
        Up = Upw.*double(image); 
        Up = parallelMove(Up,[mpc(k,1)-i,mpc(k,2)-j]); 
        oUp = oUpw.*double(image); 
        SV = (sum(sum(sum((Up-oUp).^2))))/sum(sum(Upatch(:,:,k))); 
        
        if SV<minSV 
            minSV = SV; 
            Rpatch(:,:,k) = Dpw;
            ppc(k,:) = [i,j]; 
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
    