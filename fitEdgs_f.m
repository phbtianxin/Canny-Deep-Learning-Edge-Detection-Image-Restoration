function [FittedEdgs,FittedImage] = fitEdgs_f(image,edgPm,BObj,M,N)

% [edgPm,w,edgs,M,N,Q] = getEdge_f('planet01.png');
% BObj = optimizeByGA_f(edgPm,M,N);
% image = imread('planet01.png');
% [FittedEdgs,FittedImage] = fitEdgs_f(image,edgPm,BObj,M,N);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FittedEdgs = zeros(M,N);
image = double(image);
for i=1:BObj(end) 
    
    xy1 = [edgPm(BObj(2*i-1)).Curve(:,2),edgPm(BObj(2*i-1)).Curve(:,1)]; 
    xy2 = [edgPm(BObj(2*i)).Curve(:,2),edgPm(BObj(2*i)).Curve(:,1)]; 
    p1 = [edgPm(BObj(2*i-1)).Beginning(2),edgPm(BObj(2*i-1)).Beginning(1)]; 
    p2 = [edgPm(BObj(2*i)).Beginning(2),edgPm(BObj(2*i)).Beginning(1)]; 
    
    figure(i)
    subplot(2,2,1)
    plot(xy1(:,1),xy1(:,2),'.',xy2(:,1),xy2(:,2),'.')
    title('Original discrete curve','FontSize',10,'FontAngle','italic');
    axis equal 
    
    v = [p2(1)-p1(1),p2(2)-p1(2)];
    angle = atan(v(2)/v(1)); 
    newXY1 = xy1*[cos(angle),-sin(angle);sin(angle),cos(angle)]; 
    newXY2 = xy2*[cos(angle),-sin(angle);sin(angle),cos(angle)]; 
    newp1 = newXY1(1,:); 
    newp2 = newXY2(1,:);
    
    mean5XY1 = (newXY1(5:end,1)+newXY1(4:end-1,1)+newXY1(3:end-2,1)+newXY1(2:end-3,1)+newXY1(1:end-4,1))/5; 
    dM5XY1 = mean5XY1(1:end-1) - mean5XY1(2:end); 
    position1 = find(fix(10*dM5XY1)==0);
    if length(position1)~=0
        newXY1 = newXY1(1:position1(1)+4,:);
    end
    
    mean5XY2 = (newXY2(5:end,1)+newXY2(4:end-1,1)+newXY2(3:end-2,1)+newXY2(2:end-3,1)+newXY2(1:end-4,1))/5; 
    dM5XY2 = mean5XY2(1:end-1) - mean5XY2(2:end); 
    position2 = find(fix(10*dM5XY2)==0); 
    if length(position2)~=0
        newXY2 = newXY2(1:position2(1)+4,:); 
    end
     
    figure(i)
    subplot(2,2,2)
    plot(newXY1(:,1),newXY1(:,2),'.',newXY2(:,1),newXY2(:,2),'.')
    title('Unitary discrete curve','FontSize',10,'FontAngle','italic');
    axis equal 
    
    f = fit([newXY1(:,1);newXY2(:,1)],[newXY1(:,2);newXY2(:,2)],'poly5'); 
    X = newp1(1):(newp2(1)-newp1(1))/(10*abs(newp2(1)-newp1(1))):newp2(1); 
    Y = f(X); 

    figure(i)
    subplot(2,2,3)
    plot(newXY1(:,1),newXY1(:,2),'.',newXY2(:,1),newXY2(:,2),'.')
    hold on
    plot(X,Y,'r.')
    title('Fitted curve','FontSize',10,'FontAngle','italic');
    axis equal 
     
    oXY = [X',Y]*[cos(angle),sin(angle);-sin(angle),cos(angle)]; 
    oXY = round(oXY); 
     
    figure(i)
    subplot(2,2,4)
    plot(xy1(:,1),xy1(:,2),'.',xy2(:,1),xy2(:,2),'.')
    hold on
    plot(oXY(:,1),oXY(:,2),'r.')
    title('Recovered fitted curve','FontSize',10,'FontAngle','italic');
    axis equal 
    
    FittedEdgs(sub2ind(size(FittedEdgs),oXY(:,2),oXY(:,1))) = 1; 
    
    reducedXY = [oXY(find((sum((abs(oXY(1:end-1,:)-oXY(2:end,:)))'))'~=0),:);oXY(end,:)]; 

    for l = 1:size(image,3)
        tstep = (image(reducedXY(1,2),reducedXY(1,1),l) - image(reducedXY(end,2),reducedXY(end,1),l))...
            /(size(reducedXY,1)-1); 
        for j = 2:size(reducedXY,1)-1
            image(reducedXY(j,2),reducedXY(j,1),l) = image(reducedXY(j-1,2),reducedXY(j-1,1),l) - tstep; 
        end
    end
end
FittedEdgs = FittedEdgs==1; 
FittedImage = uint8(image); 

imshow(FittedImage)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     