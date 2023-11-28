function [edgPm,w,edgs,M,N,Q] = getEdge_f(image)
%��ȡ����������ͼ���������������нӴ��ĸ������������м���Ӧ�������

%���У�
% edgPm --- һ���ṹ�����У��м������������Ӵ��������߶˵���м��������Ľṹ�壩��ÿ���ṹ��������������
% 1������������Ӵ���һ�������ߵĶ˵�����
% 2���ö˵������������ߵ���������
% 3���ò������������������ƽ������
% 4���ò�����������������ĶԱȶȣ����
% image --- һ���ַ���������һ��ͼ��
% w --- ͼ�������������Ĥ
% edgs --- ������������ߵ���Ĥ
% [M,N,Q] --- ͼ��ĳߴ�

%ʾ����
% [edgPm,w,edgs,M,N,Q] = getEdge_f('planet01.png');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%image = imread(image); %��ȡͼ�����
[M,N,Q] = size(image); %���ͼ��ĳߴ�
if Q==3
    image = rgb2gray(image); %���ͼ��ĻҶ�ֵ����
end

w = image==0; %���ͼ�������������Ĥ

%�������ͼ����������2ֵ����
g = edge(image, 'canny', [0.4, 0.5], 1.4); 
% 'canny'���������е�����ֵԽ���ҳ�������Խǿ�����һ��ֵԽ���ҳ�������Խƽ��
imshow(g);
%����ǰ����������Ĥ
w_N = [w(2:M,:);zeros(1,N)]==1; %�ϴ�һ��
w_S = [zeros(1,N);w(1:M-1,:)]==1; %�´�һ��
w_W = [w(:,2:N),zeros(M,1)]==1; %���һ��
w_E = [zeros(M,1),w(:,1:N-1)]==1; %�Ҵ�һ��
w_NW = [w_N(:,2:N),zeros(M,1)]==1; %�ϴ�һ�в����һ��
w_NE = [zeros(M,1),w_N(:,1:N-1)]==1; %�ϴ�һ�в��Ҵ�һ��
w_SW = [w_S(:,2:N),zeros(M,1)]==1; %�´�һ�в����һ��
w_SE = [zeros(M,1),w_S(:,1:N-1)]==1; %�´�һ�в��Ҵ�һ��
bw = w_N|w_S|w_W|w_E|w_NW|w_NE|w_SW|w_SE|w; %����������Ĥ��һȦ��
newg = (g-bw)==1; %ȥ����������Ĥ������

%����ǰ����������Ĥ
w_N = [bw(2:M,:);zeros(1,N)]==1; %�ϴ�һ��
w_S = [zeros(1,N);bw(1:M-1,:)]==1; %�´�һ��
w_W = [bw(:,2:N),zeros(M,1)]==1; %���һ��
w_E = [zeros(M,1),bw(:,1:N-1)]==1; %�Ҵ�һ��
w_NW = [w_N(:,2:N),zeros(M,1)]==1; %�ϴ�һ�в����һ��
w_NE = [zeros(M,1),w_N(:,1:N-1)]==1; %�ϴ�һ�в��Ҵ�һ��
w_SW = [w_S(:,2:N),zeros(M,1)]==1; %�´�һ�в����һ��
w_SE = [zeros(M,1),w_S(:,1:N-1)]==1; %�´�һ�в��Ҵ�һ��
bbw = w_N|w_S|w_W|w_E|w_NW|w_NE|w_SW|w_SE|bw; %����������Ĥ����һȦ��

%����ǰ����������Ĥ
w_N = [bbw(2:M,:);zeros(1,N)]==1; %�ϴ�һ��
w_S = [zeros(1,N);bbw(1:M-1,:)]==1; %�´�һ��
w_W = [bbw(:,2:N),zeros(M,1)]==1; %���һ��
w_E = [zeros(M,1),bbw(:,1:N-1)]==1; %�Ҵ�һ��
w_NW = [w_N(:,2:N),zeros(M,1)]==1; %�ϴ�һ�в����һ��
w_NE = [zeros(M,1),w_N(:,1:N-1)]==1; %�ϴ�һ�в��Ҵ�һ��
w_SW = [w_S(:,2:N),zeros(M,1)]==1; %�´�һ�в����һ��
w_SE = [zeros(M,1),w_S(:,1:N-1)]==1; %�´�һ�в��Ҵ�һ��
bbbw = w_N|w_S|w_W|w_E|w_NW|w_NE|w_SW|w_SE|bbw; %����������Ĥ����һȦ��

%ѡ�������������нӴ���������
tbind = find(bw==1); %�ҵ��������е�һ�����ص�����
tL = bwlabel(bbw|newg,8); %������������Ӵ��������߼��������ϲ���һ�����ӷ���
tL = tL==tL(tbind(1)); %�������������ӷ�����ֵ���1������ȫ��0���Ӷ��õ�һ��ֻ����������������������Ӵ������ߵ���Ĥ
finalg = newg&tL; %��������Ĥ�ڳ�����������Ӵ���������
tL = bwlabel(finalg,8); %��ʾ�����������
for i=1:max(max(tL))
    [tI,tJ] = find(tL==i);
    edgeCurve{i,1} = [tI,tJ];
end %����ÿ�������ߵ���������

edgs = zeros(M,N); %��ʼ�������������Ĥ
for i=1:length(edgeCurve)
    edgs(sub2ind([M,N],edgeCurve{i,1}(:,1),edgeCurve{i,1}(:,2))) = 1;
end %��������������ߵ�������1
edgs = edgs==1; %2ֵ��

%%��������ȡ��������ͼ
%imview(edgs)
%imview(w)
imshow(edgs)

%������������ߵĸ������
d = 5; %�������ȡ��Ŀ��
k = 1;
AllC = [];
for i=1:length(edgeCurve)
    finalg = zeros(M,N); %��ʼ����i����������Ĥ
    finalg(sub2ind([M,N], edgeCurve{i,1}(:,1), edgeCurve{i,1}(:,2)))=1; %�����i��������
    finalg = finalg==1; %2ֵ��
    
    n = 1;
    for j=1:size(edgeCurve{i,1},1) %����ĳ�������ϵ���������
        if sum(sum(finalg(edgeCurve{i,1}(j,1)-1:edgeCurve{i,1}(j,1)+1,edgeCurve{i,1}(j,2)-1:edgeCurve{i,1}(j,2)+1)))<3
            Cind(n,:) = edgeCurve{i,1}(j,:); %������ؼ�Ϊ�˵�����
            %��ĳһ������Χ��֮����������ֻ��һ��
            n = n+1;
        elseif sum(sum(finalg(edgeCurve{i,1}(j,1)-1:edgeCurve{i,1}(j,1)+1,edgeCurve{i,1}(j,2)-1:edgeCurve{i,1}(j,2)+1)))==3
             if (sum(sum(finalg(edgeCurve{i,1}(j,1)-1:edgeCurve{i,1}(j,1),edgeCurve{i,1}(j,2)-1:edgeCurve{i,1}(j,2)+1)))==1)||(sum(sum(finalg(edgeCurve{i,1}(j,1):edgeCurve{i,1}(j,1)+1,edgeCurve{i,1}(j,2)-1:edgeCurve{i,1}(j,2)+1)))==1)
                  %delete middle point of the edgeCurve 
                  if (edgeCurve{i,1}(j,2)-1==1)&&(sum(finalg(edgeCurve{i,1}(j,1),edgeCurve{i,1}(j,2):edgeCurve{i,1}(j,2)+2))<2)
                       Cind(n,:) = edgeCurve{i,1}(j,:); %������ؼ�Ϊ�˵�����
                       n = n+1;
                  elseif (edgeCurve{i,1}(j,2)+1==N)&&(sum(finalg(edgeCurve{i,1}(j,1),edgeCurve{i,1}(j,2)-2:edgeCurve{i,1}(j,2)))<2)
                       Cind(n,:) = edgeCurve{i,1}(j,:); %������ؼ�Ϊ�˵�����
                       n = n+1;
                  elseif (edgeCurve{i,1}(j,2)-1>1)&&(edgeCurve{i,1}(j,2)+1<N)
                      if (sum(finalg(edgeCurve{i,1}(j,1),edgeCurve{i,1}(j,2)-2:edgeCurve{i,1}(j,2)))<2)&&(sum(finalg(edgeCurve{i,1}(j,1),edgeCurve{i,1}(j,2):edgeCurve{i,1}(j,2)+2))<2)
                          Cind(n,:) = edgeCurve{i,1}(j,:); %������ؼ�Ϊ�˵�����
                          n = n+1;
                      else
                          finalg(edgeCurve{i,1}(j,1),edgeCurve{i,1}(j,2))=0;
                      end  
                  end                  
             end
             if (sum(sum(finalg(edgeCurve{i,1}(j,1)-1:edgeCurve{i,1}(j,1)+1,edgeCurve{i,1}(j,2):edgeCurve{i,1}(j,2)+1)))==1)||(sum(sum(finalg(edgeCurve{i,1}(j,1)-1:edgeCurve{i,1}(j,1)+1,edgeCurve{i,1}(j,2)-1:edgeCurve{i,1}(j,2))))==1)
                 %delete middle point of the edgeCurve 
                 if (edgeCurve{i,1}(j,1)-1==1)&&(sum(finalg(edgeCurve{i,1}(j,1):edgeCurve{i,1}(j,1)+2,edgeCurve{i,1}(j,2)))<2)
                     Cind(n,:) = edgeCurve{i,1}(j,:); %������ؼ�Ϊ�˵�����
                     n = n+1;
                 elseif (edgeCurve{i,1}(j,1)+1==M)&&(sum(finalg(edgeCurve{i,1}(j,1)-2:edgeCurve{i,1}(j,1),edgeCurve{i,1}(j,2)))<2)
                     Cind(n,:) = edgeCurve{i,1}(j,:); %������ؼ�Ϊ�˵�����
                     n = n+1;                     
                 elseif (edgeCurve{i,1}(j,1)-1>1)&&(edgeCurve{i,1}(j,1)+1<M)
                     if (sum(finalg(edgeCurve{i,1}(j,1)-2:edgeCurve{i,1}(j,1),edgeCurve{i,1}(j,2)))<2)&&(sum(finalg(edgeCurve{i,1}(j,1):edgeCurve{i,1}(j,1)+2,edgeCurve{i,1}(j,2)))<2)
                         Cind(n,:) = edgeCurve{i,1}(j,:); %������ؼ�Ϊ�˵�����
                         n = n+1;                                         
                     else
                         finalg(edgeCurve{i,1}(j,1),edgeCurve{i,1}(j,2))=0;
                     end
                 end
             end
        end
    end%��õ�i�����������˵�����

    %%Ϊ���Ƹ���������ߵĶ˵�
    %AllC = [AllC;Cind];
    
    if bbbw(Cind(1,1),Cind(1,2)) %����˵�1����������Ӵ�
        pw = zeros(M,N); %��ʼ���˵��������ȡ����Ĥ
        
        tw = pw; %Ϊ����������϶˵��ǰһ�������ʼ��һ����ʱ��Ĥ
        tw(Cind(1,1)-1:Cind(1,1)+1,Cind(1,2)-1:Cind(1,2)+1) = 1; %�ڶ˵㴦��һȦ��
        tw(Cind(1,1),Cind(1,2)) = 0; %��ȥ�˵�
        tw = tw==1; %2ֵ��
        tw = tw&finalg; %�ڳ��������϶˵��ǰһ��
        [R,C] = find(tw==1); %�õ������϶˵��ǰһ�������
        %Two point method
        if size([R,C],1)>1
            R=R(1);
            C=C(1);
        end  
        
        anticR = Cind(1,1) + round(((R-Cind(1,1))-(C-Cind(1,2)))/2); %�õ��Զ˵�Ϊ������ʱ��ת����һ���������
        anticC = Cind(1,2) + round(((R-Cind(1,1))+(C-Cind(1,2)))/2); %�õ��Զ˵�Ϊ������ʱ��ת����һ���������
        
        pw(Cind(1,1)-d:Cind(1,1)+d,Cind(1,2)-d:Cind(1,2)+d) = 1; %������������������֮��ľ��붼�����������ȡ��Ŀ�ȣ�����������Σ���
        pw = pw==1; %�õ��������߶˵�Ϊ���ĵ�һ�����δ�����Ĥ
        
        region = pw&~(bbw|finalg); %������δ����г������߼�һ������������Ĳ���
       
        L = bwlabel(region,4); %�õ����δ������������߼�һ�������������ֿ����������ӷ���
        pwRight = L==L(anticR,anticC); %��ô�����������������ȥ���Ҳಿ�ֵ���ȡ��Ĥ
        pwLeft = region&(~pwRight); %��ô�����������������ȥ��ಿ�ֵ���ȡ��Ĥ
        
        %%���ƴ�i������������������ȥ�������ಿ�ֵ���ȡ��Ĥ
        %imview(pwRight);
        %imview(pwLeft);
        
        edgPm(k).Beginning = [Cind(1,1),Cind(1,2)]; %����i�����������������ĽӴ��������¼������ṹ��Beginning�������
        
        %����������Ը������ߵ���������������
        finalg(edgPm(k).Beginning(1),edgPm(k).Beginning(2)) = 0; %ȥ���˵�
        clear lineSub

        [lineSub(:,1),lineSub(:,2)] = find(finalg==1); %�ҵ��������е�����
        newCurve = edgPm(k).Beginning; %���˵��������������еĵ�һ��
        while size(lineSub,1)~=0
            distance = abs(ones(size(lineSub,1),1)*newCurve(end,:)-lineSub); %����ǰһ�������¸���ľ���
            sumd = distance(:,1).^2 + distance(:,2).^2; %�õ����������ֵ��ƽ����
            if find(sumd<=1)
                newCurve = [newCurve;lineSub(find(sumd<=1),:)]; %���ҵ�����һ��������ӵ���������
            end
            if find(sumd>1&sumd<3)
                newCurve = [newCurve;lineSub(find(sumd>1&sumd<3),:)]; %���ҵ�����һ��������ӵ���������            
            end   
            lineSub(find(sumd<3),:) = []; %�ھ���������ȥ�õ�����
        end
        
        edgPm(k).Curve = newCurve; %����i�����������м�¼������ṹ��Curve�������
        edgPm(k).MeanColor = [sum(sum(pwLeft.*double(image)))/sum(sum(pwLeft)),sum(sum(pwRight.*double(image)))/sum(sum(pwRight))]; %%
        %���ӵ�i������������������ȥ����������ͼ���ƽ�����ȼ�¼������ṹ��MeanColor�������
        edgPm(k).variance = [sqrt((sum(sum((pwLeft.*(double(image)-edgPm(k).MeanColor(1,1))).^2)))/sum(sum(pwLeft))),...
            sqrt((sum(sum((pwRight.*(double(image)-edgPm(k).MeanColor(1,2))).^2)))/sum(sum(pwRight)))];
        k = k+1;
    end
    
    if bbbw(Cind(2,1),Cind(2,2)) %����˵�2����������Ӵ�
        pw = zeros(M,N); %��ʼ���˵��������ȡ����Ĥ
        
        tw = pw; %Ϊ����������϶˵��ǰһ�������ʼ��һ����ʱ��Ĥ
        tw(Cind(2,1)-1:Cind(2,1)+1,Cind(2,2)-1:Cind(2,2)+1) = 1; %�ڶ˵㴦��һȦ��
        tw(Cind(2,1),Cind(2,2)) = 0; %��ȥ�˵�
        tw = tw==1; %2ֵ��
        tw = tw&finalg; %�ڳ��������϶˵��ǰһ��
        [R,C] = find(tw==1); %�õ������϶˵��ǰһ�������
        %Two point method
        if size([R,C],1)>1
            R=R(1);
            C=C(1);
        end               
        anticR = Cind(2,1) + round(((R-Cind(2,1))-(C-Cind(2,2)))/2); %�õ��Զ˵�Ϊ������ʱ��ת����һ���������
        anticC = Cind(2,2) + round(((R-Cind(2,1))+(C-Cind(2,2)))/2); %�õ��Զ˵�Ϊ������ʱ��ת����һ���������
        
        pw(Cind(2,1)-d:Cind(2,1)+d,Cind(2,2)-d:Cind(2,2)+d) = 1; %������������������֮��ľ��붼�����������ȡ��Ŀ��
        pw = pw==1; %�õ��������߶˵�Ϊ���ĵ�һ�����δ�����Ĥ
        
        region = pw&~(bbw|finalg); %������δ����г������߼�һ������������Ĳ���
        
        L = bwlabel(region,4); %�õ����δ������������߼�һ�������������ֿ����������ӷ���
        pwRight = L==L(anticR,anticC); %��ô�����������������ȥ���Ҳಿ�ֵ���ȡ��Ĥ
        pwLeft = region&(~pwRight); %��ô�����������������ȥ��ಿ�ֵ���ȡ��Ĥ
        
        %%���ƴ�i������������������ȥ�������ಿ�ֵ���ȡ��Ĥ
        %imview(pwRight);
        %imview(pwLeft);
        
        edgPm(k).Beginning = [Cind(2,1),Cind(2,2)]; %����i�����������������ĽӴ��������¼������ṹ��Beginning�������
        
        %����������Ը������ߵ���������������
        finalg(edgPm(k).Beginning(1),edgPm(k).Beginning(2)) = 0; %ȥ���˵�
        clear lineSub
        [lineSub(:,1),lineSub(:,2)] = find(finalg==1); %�ҵ��������е�����
        newCurve = edgPm(k).Beginning; %���˵��������������еĵ�һ��
        while size(lineSub,1)~=0
            distance = abs(ones(size(lineSub,1),1)*newCurve(end,:)-lineSub); %����ǰһ�������¸���ľ���
            sumd = distance(:,1).^2 + distance(:,2).^2; %�õ����������ֵ��ƽ����
            %write oder rom small to large
            if find(sumd<=1)
                newCurve = [newCurve;lineSub(find(sumd<=1),:)]; %���ҵ�����һ��������ӵ���������
            end
            if find(sumd>1&sumd<3)
                newCurve = [newCurve;lineSub(find(sumd>1&sumd<3),:)]; %���ҵ�����һ��������ӵ���������            
            end            
            lineSub(find(sumd<3),:) = []; %�ھ���������ȥ�õ�����
        end
        
        edgPm(k).Curve = newCurve; %����i�����������м�¼������ṹ��Curve�������
        edgPm(k).MeanColor = [sum(sum(pwLeft.*double(image)))/sum(sum(pwLeft)),sum(sum(pwRight.*double(image)))/sum(sum(pwRight))]; %%
        %���ӵ�i������������������ȥ����������ͼ���ƽ�����ȼ�¼������ṹ��MeanColor�������
        edgPm(k).variance = [sqrt((sum(sum((pwLeft.*(double(image)-edgPm(k).MeanColor(1,1))).^2)))/sum(sum(pwLeft))),...
            sqrt((sum(sum((pwRight.*(double(image)-edgPm(k).MeanColor(1,2))).^2)))/sum(sum(pwRight)))];
        k = k+1;
    end
end

%%���Ƹ���������ߵĶ˵�
%drawp = zeros(M,N); %��ʼ��
%drawp(sub2ind([M,N], AllC(:,1), AllC(:,2))) = 1; %���˵��1
%drawp = drawp==1; %2ֵ��
%imview(drawp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


