function [edgPm,w,edgs,M,N,Q] = getEdge_f(image)
% [edgPm,w,edgs,M,N,Q] = getEdge_f('planet01.png');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[M,N,Q] = size(image); 
if Q==3
    image = rgb2gray(image); 
end

w = image==0; 
g = edge(image, 'canny', [0.4, 0.5], 1.4); 
imshow(g);

w_N = [w(2:M,:);zeros(1,N)]==1; 
w_S = [zeros(1,N);w(1:M-1,:)]==1; 
w_W = [w(:,2:N),zeros(M,1)]==1; 
w_E = [zeros(M,1),w(:,1:N-1)]==1; 
w_NW = [w_N(:,2:N),zeros(M,1)]==1; 
w_NE = [zeros(M,1),w_N(:,1:N-1)]==1; 
w_SW = [w_S(:,2:N),zeros(M,1)]==1; 
w_SE = [zeros(M,1),w_S(:,1:N-1)]==1; 
bw = w_N|w_S|w_W|w_E|w_NW|w_NE|w_SW|w_SE|w; 
newg = (g-bw)==1; 

w_N = [bw(2:M,:);zeros(1,N)]==1; 
w_S = [zeros(1,N);bw(1:M-1,:)]==1; 
w_W = [bw(:,2:N),zeros(M,1)]==1; 
w_E = [zeros(M,1),bw(:,1:N-1)]==1; 
w_NW = [w_N(:,2:N),zeros(M,1)]==1; 
w_NE = [zeros(M,1),w_N(:,1:N-1)]==1; 
w_SW = [w_S(:,2:N),zeros(M,1)]==1; 
w_SE = [zeros(M,1),w_S(:,1:N-1)]==1; 
bbw = w_N|w_S|w_W|w_E|w_NW|w_NE|w_SW|w_SE|bw; 

w_N = [bbw(2:M,:);zeros(1,N)]==1; 
w_S = [zeros(1,N);bbw(1:M-1,:)]==1; 
w_W = [bbw(:,2:N),zeros(M,1)]==1; 
w_E = [zeros(M,1),bbw(:,1:N-1)]==1; 
w_NW = [w_N(:,2:N),zeros(M,1)]==1; 
w_NE = [zeros(M,1),w_N(:,1:N-1)]==1; 
w_SW = [w_S(:,2:N),zeros(M,1)]==1; 
w_SE = [zeros(M,1),w_S(:,1:N-1)]==1; 
bbbw = w_N|w_S|w_W|w_E|w_NW|w_NE|w_SW|w_SE|bbw; 

tbind = find(bw==1); 
tL = bwlabel(bbw|newg,8); 
tL = tL==tL(tbind(1)); 
finalg = newg&tL; 
tL = bwlabel(finalg,8); 
for i=1:max(max(tL))
    [tI,tJ] = find(tL==i);
    edgeCurve{i,1} = [tI,tJ];
end 

edgs = zeros(M,N); 
for i=1:length(edgeCurve)
    edgs(sub2ind([M,N],edgeCurve{i,1}(:,1),edgeCurve{i,1}(:,2))) = 1;
end 
edgs = edgs==1; 

imshow(edgs)

d = 5; 
k = 1;
AllC = [];
for i=1:length(edgeCurve)
    finalg = zeros(M,N); 
    finalg(sub2ind([M,N], edgeCurve{i,1}(:,1), edgeCurve{i,1}(:,2)))=1; 
    finalg = finalg==1; 
    
    n = 1;
    for j=1:size(edgeCurve{i,1},1) 
        if sum(sum(finalg(edgeCurve{i,1}(j,1)-1:edgeCurve{i,1}(j,1)+1,edgeCurve{i,1}(j,2)-1:edgeCurve{i,1}(j,2)+1)))<3
            Cind(n,:) = edgeCurve{i,1}(j,:); 

            n = n+1;
        elseif sum(sum(finalg(edgeCurve{i,1}(j,1)-1:edgeCurve{i,1}(j,1)+1,edgeCurve{i,1}(j,2)-1:edgeCurve{i,1}(j,2)+1)))==3
             if (sum(sum(finalg(edgeCurve{i,1}(j,1)-1:edgeCurve{i,1}(j,1),edgeCurve{i,1}(j,2)-1:edgeCurve{i,1}(j,2)+1)))==1)||(sum(sum(finalg(edgeCurve{i,1}(j,1):edgeCurve{i,1}(j,1)+1,edgeCurve{i,1}(j,2)-1:edgeCurve{i,1}(j,2)+1)))==1)
 
                  if (edgeCurve{i,1}(j,2)-1==1)&&(sum(finalg(edgeCurve{i,1}(j,1),edgeCurve{i,1}(j,2):edgeCurve{i,1}(j,2)+2))<2)
                       Cind(n,:) = edgeCurve{i,1}(j,:); 
                       n = n+1;
                  elseif (edgeCurve{i,1}(j,2)+1==N)&&(sum(finalg(edgeCurve{i,1}(j,1),edgeCurve{i,1}(j,2)-2:edgeCurve{i,1}(j,2)))<2)
                       Cind(n,:) = edgeCurve{i,1}(j,:); 
                       n = n+1;
                  elseif (edgeCurve{i,1}(j,2)-1>1)&&(edgeCurve{i,1}(j,2)+1<N)
                      if (sum(finalg(edgeCurve{i,1}(j,1),edgeCurve{i,1}(j,2)-2:edgeCurve{i,1}(j,2)))<2)&&(sum(finalg(edgeCurve{i,1}(j,1),edgeCurve{i,1}(j,2):edgeCurve{i,1}(j,2)+2))<2)
                          Cind(n,:) = edgeCurve{i,1}(j,:); 
                          n = n+1;
                      else
                          finalg(edgeCurve{i,1}(j,1),edgeCurve{i,1}(j,2))=0;
                      end  
                  end                  
             end
             if (sum(sum(finalg(edgeCurve{i,1}(j,1)-1:edgeCurve{i,1}(j,1)+1,edgeCurve{i,1}(j,2):edgeCurve{i,1}(j,2)+1)))==1)||(sum(sum(finalg(edgeCurve{i,1}(j,1)-1:edgeCurve{i,1}(j,1)+1,edgeCurve{i,1}(j,2)-1:edgeCurve{i,1}(j,2))))==1)
 
                 if (edgeCurve{i,1}(j,1)-1==1)&&(sum(finalg(edgeCurve{i,1}(j,1):edgeCurve{i,1}(j,1)+2,edgeCurve{i,1}(j,2)))<2)
                     Cind(n,:) = edgeCurve{i,1}(j,:);
                     n = n+1;
                 elseif (edgeCurve{i,1}(j,1)+1==M)&&(sum(finalg(edgeCurve{i,1}(j,1)-2:edgeCurve{i,1}(j,1),edgeCurve{i,1}(j,2)))<2)
                     Cind(n,:) = edgeCurve{i,1}(j,:); 
                     n = n+1;                     
                 elseif (edgeCurve{i,1}(j,1)-1>1)&&(edgeCurve{i,1}(j,1)+1<M)
                     if (sum(finalg(edgeCurve{i,1}(j,1)-2:edgeCurve{i,1}(j,1),edgeCurve{i,1}(j,2)))<2)&&(sum(finalg(edgeCurve{i,1}(j,1):edgeCurve{i,1}(j,1)+2,edgeCurve{i,1}(j,2)))<2)
                         Cind(n,:) = edgeCurve{i,1}(j,:); 
                         n = n+1;                                         
                     else
                         finalg(edgeCurve{i,1}(j,1),edgeCurve{i,1}(j,2))=0;
                     end
                 end
             end
        end
    end
    
    if bbbw(Cind(1,1),Cind(1,2)) 
        pw = zeros(M,N); 
        
        tw = pw; 
        tw(Cind(1,1)-1:Cind(1,1)+1,Cind(1,2)-1:Cind(1,2)+1) = 1;
        tw(Cind(1,1),Cind(1,2)) = 0; 
        tw = tw==1; 
        tw = tw&finalg; 
        [R,C] = find(tw==1); 

        if size([R,C],1)>1
            R=R(1);
            C=C(1);
        end  
        
        anticR = Cind(1,1) + round(((R-Cind(1,1))-(C-Cind(1,2)))/2); 
        anticC = Cind(1,2) + round(((R-Cind(1,1))+(C-Cind(1,2)))/2); 
        
        pw(Cind(1,1)-d:Cind(1,1)+d,Cind(1,2)-d:Cind(1,2)+d) = 1; 
        pw = pw==1; 
        
        region = pw&~(bbw|finalg); 
       
        L = bwlabel(region,4); 
        pwRight = L==L(anticR,anticC); 
        pwLeft = region&(~pwRight); 
        
        edgPm(k).Beginning = [Cind(1,1),Cind(1,2)]; 
        
        finalg(edgPm(k).Beginning(1),edgPm(k).Beginning(2)) = 0; 
        clear lineSub

        [lineSub(:,1),lineSub(:,2)] = find(finalg==1); 
        newCurve = edgPm(k).Beginning; 
        while size(lineSub,1)~=0
            distance = abs(ones(size(lineSub,1),1)*newCurve(end,:)-lineSub); 
            sumd = distance(:,1).^2 + distance(:,2).^2; 
            if find(sumd<=1)
                newCurve = [newCurve;lineSub(find(sumd<=1),:)]; 
            end
            if find(sumd>1&sumd<3)
                newCurve = [newCurve;lineSub(find(sumd>1&sumd<3),:)];            
            end   
            lineSub(find(sumd<3),:) = []; 
        end
        
        edgPm(k).Curve = newCurve; 
        edgPm(k).MeanColor = [sum(sum(pwLeft.*double(image)))/sum(sum(pwLeft)),sum(sum(pwRight.*double(image)))/sum(sum(pwRight))]; 
        edgPm(k).variance = [sqrt((sum(sum((pwLeft.*(double(image)-edgPm(k).MeanColor(1,1))).^2)))/sum(sum(pwLeft))),...
            sqrt((sum(sum((pwRight.*(double(image)-edgPm(k).MeanColor(1,2))).^2)))/sum(sum(pwRight)))];
        k = k+1;
    end
    
    if bbbw(Cind(2,1),Cind(2,2)) 
        pw = zeros(M,N); 
        
        tw = pw; 
        tw(Cind(2,1)-1:Cind(2,1)+1,Cind(2,2)-1:Cind(2,2)+1) = 1; 
        tw(Cind(2,1),Cind(2,2)) = 0; 
        tw = tw==1; 
        tw = tw&finalg; 
        [R,C] = find(tw==1); 
        if size([R,C],1)>1
            R=R(1);
            C=C(1);
        end               
        anticR = Cind(2,1) + round(((R-Cind(2,1))-(C-Cind(2,2)))/2); 
        anticC = Cind(2,2) + round(((R-Cind(2,1))+(C-Cind(2,2)))/2); 
        
        pw(Cind(2,1)-d:Cind(2,1)+d,Cind(2,2)-d:Cind(2,2)+d) = 1; 
        pw = pw==1; 
        
        region = pw&~(bbw|finalg); 
        
        L = bwlabel(region,4); 
        pwRight = L==L(anticR,anticC); 
        pwLeft = region&(~pwRight); 
        
        edgPm(k).Beginning = [Cind(2,1),Cind(2,2)]; 
        
        finalg(edgPm(k).Beginning(1),edgPm(k).Beginning(2)) = 0; 
        clear lineSub
        [lineSub(:,1),lineSub(:,2)] = find(finalg==1); 
        newCurve = edgPm(k).Beginning; 
        while size(lineSub,1)~=0
            distance = abs(ones(size(lineSub,1),1)*newCurve(end,:)-lineSub); 
            sumd = distance(:,1).^2 + distance(:,2).^2; 
            if find(sumd<=1)
                newCurve = [newCurve;lineSub(find(sumd<=1),:)]; 
            end
            if find(sumd>1&sumd<3)
                newCurve = [newCurve;lineSub(find(sumd>1&sumd<3),:)];         
            end            
            lineSub(find(sumd<3),:) = []; 
        end
        
        edgPm(k).Curve = newCurve; 
        edgPm(k).MeanColor = [sum(sum(pwLeft.*double(image)))/sum(sum(pwLeft)),sum(sum(pwRight.*double(image)))/sum(sum(pwRight))]; %%
        edgPm(k).variance = [sqrt((sum(sum((pwLeft.*(double(image)-edgPm(k).MeanColor(1,1))).^2)))/sum(sum(pwLeft))),...
            sqrt((sum(sum((pwRight.*(double(image)-edgPm(k).MeanColor(1,2))).^2)))/sum(sum(pwRight)))];
        k = k+1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


