function RImage = GAmain_f(Image)
% RImage = GAmain_f('planet01.png');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

image=Image;
imshow(image);

[edgPm,w,edgs,M,N,Q] = getEdge_f(Image);

BObj = optimizeByGA_f(edgPm,M,N);

[FittedEdgs,FittedImage] = fitEdgs_f(image,edgPm,BObj,M,N);

RImage = double(FittedImage); 

Record=[(sum(sum(w&~FittedEdgs))),1];
while sum(sum(w&~FittedEdgs))~=0 
    
    [mpc,Upatch,Dpatch,tp] = GAfindMPC_f(w,edgs,FittedEdgs);
    
    [ppc,Rpatch] = GAfindPP_f(RImage,w,edgs,FittedEdgs,mpc,Upatch,Dpatch,tp);
    
    for i=1:size(ppc,1)
        Rpw = Rpatch(:,:,i); 
        
        if Q==3 
            Rpw = cat(3,Rpw,cat(3,Rpw,Rpw)); 
        end
        
        Rp = Rpw.*RImage; 
        Rp = parallelMove(Rp,[mpc(i,1)-ppc(i,1),mpc(i,2)-ppc(i,2)]); 
        RImage = RImage + Rp; 
        w = w&(~(Rp(:,:,1)>0));
    end

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
RImage = uint8(RImage); 

imshow(RImage);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        