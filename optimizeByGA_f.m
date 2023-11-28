function BObj = optimizeByGA_f(edgPm,M,N)

% [edgPm,w,edgs,M,N,Q] = getEdge_f('planet01.png');
% BObj = optimizeByGA_f(edgPm,M,N);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nind = 20; 
MaxGen = 100; 
T = 1.2; 
Lind = size(edgPm,2); 
BaseV = [Lind:-1:1,floor(Lind/2)]; 
Chrom = crtbp(Nind,BaseV); 
Pm = 0.7/Lind; 
XOVR = 0.3; 

traceM = []; 
traceB = []; 

difference = zeros(Lind); 
    
figurek = 1;
Tdifference = 0*difference; 
for i=1:Lind-1
    for j=i+1:Lind 
        
        xy1 = [edgPm(i).Curve(:,2),edgPm(i).Curve(:,1)]; 
        xy2 = [edgPm(j).Curve(:,2),edgPm(j).Curve(:,1)]; 
        p1 = [edgPm(i).Beginning(2),edgPm(i).Beginning(1)]; 
        p2 = [edgPm(j).Beginning(2),edgPm(j).Beginning(1)]; 

        v = [p2(1)-p1(1),p2(2)-p1(2)]; 
        angle = atan(v(2)/v(1)); 
        newXY1 = xy1*[cos(angle),-sin(angle);sin(angle),cos(angle)]; 
        newXY2 = xy2*[cos(angle),-sin(angle);sin(angle),cos(angle)]; 
        newp1 = newXY1(1,:);
        newp2 = newXY2(1,:); 

        mean5XY1 = (newXY1(5:end,1)+newXY1(4:end-1,1)+newXY1(3:end-2,1)+newXY1(2:end-3,1)+newXY1(1:end-4,1))/5; 
        dM5XY1 = mean5XY1(1:end-1) - mean5XY1(2:end); 
        position1 = find(round(dM5XY1)==0); 
        if length(position1)~=0
            newXY1 = newXY1(1:position1(1)+4,:); 
        end
        
        mean5XY2 = (newXY2(5:end,1)+newXY2(4:end-1,1)+newXY2(3:end-2,1)+newXY2(2:end-3,1)+newXY2(1:end-4,1))/5; 
        dM5XY2 = mean5XY2(1:end-1) - mean5XY2(2:end); 
        position2 = find(round(dM5XY2)==0); 
        if length(position2)~=0
            newXY2 = newXY2(1:position2(1)+4,:); 
        end

        if ~(((newXY1(end,1)<newXY1(1,1))&&(newXY1(1,1)<newXY2(1,1))&&(newXY2(1,1)<newXY2(end,1)))||...
                ((newXY1(end,1)>newXY1(1,1))&&(newXY1(1,1)>newXY2(1,1))&&(newXY2(1,1)>newXY2(end,1)))) 
            Tdifference(i,j) = inf; 
            continue
        end

        f = fit([newXY1(:,1);newXY2(:,1)],[newXY1(:,2);newXY2(:,2)],'poly5');

        clear cc;
        cc(1) = f.p1; 
        cc(2) = f.p2; 
        cc(3) = f.p3; 
        cc(4) = f.p4; 
        cc(5) = f.p5;
        cc(6) = f.p6;
        x=sym('x'); 
        cf = cc*[x^5;x^4;x^3;x^2;x;1]; 
        cK = abs(diff(diff(cf)))/(1+(diff(cf))^2)^(3/2); 
        dcK = diff(cK); 
        dx = round(newp1(1)):(newp2(1)-newp1(1))/abs(newp2(1)-newp1(1)):round(newp2(1)); 
        SdcK = abs(subs(dcK,x,dx)); 
        difference(i,j) = sum(SdcK)/length(dx); 
        
        figurek = figurek + 1;
    end
end
difference = difference./(sum(sum(difference))/(sum(sum(difference>0)))); 
difference = difference + Tdifference;
    
Tdifference = 0*difference; 
for i=1:Lind
    for j=i+1:Lind 
        Tdifference(i,j) = abs(edgPm(i).MeanColor(1)-edgPm(j).MeanColor(2))+abs(edgPm(i).MeanColor(2)-edgPm(j).MeanColor(1));
    end
end
Tdifference = Tdifference./(sum(sum(Tdifference))/(sum(sum(Tdifference>0))));
difference = difference + Tdifference; 
difference = difference + difference'; 

gen = 0; 
while gen < MaxGen

    Chrom = mut(Chrom,Pm,BaseV);

    Chrom = xovsp(Chrom, XOVR);

    Chrom = Chrom+1;
    Phenotype = 0*Chrom; 
    for i=1:Nind
        tem = 1:Lind; 
        for j=1:2*Chrom(i,end) 
            Phenotype(i,j) = tem(Chrom(i,j)); 
            tem(Chrom(i,j)) = []; 
        end
    end
    Phenotype(:,end)=Chrom(:,Lind+1);
    Chrom = Chrom-1; 

    tPhenotype = Phenotype + (Phenotype==0); 
    ObjV = zeros(Nind,1); 
    for i=1:max(tPhenotype(:,Lind+1)) 
        ObjV = ObjV + difference(sub2ind(size(difference), tPhenotype(:,2*i-1), tPhenotype(:,2*i))); 
    end 
    
    ObjV = (T+ObjV)./(tPhenotype(:,Lind+1)); 
    traceM = [traceM,mean(ObjV(find(ObjV~=inf)))]; 
    traceB = [traceB,min(ObjV)];

    FitnV = ranking(ObjV); 

    NewChrIx = rws(FitnV, Nind);
    Chrom = Chrom(NewChrIx, :);

    gen = gen+1;
end

[BFV,BFI] = max(FitnV); 
BObj = Phenotype(BFI,:);

figure(figurek),
subplot(2,1,1)
plot(traceM);
title('The trace of the mean object value of population','FontSize',10,'FontAngle','italic');
subplot(2,1,2)
plot(traceB,'.');
title('The trace of the best object value of population','FontSize',10,'FontAngle','italic');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

