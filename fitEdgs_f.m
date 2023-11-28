function [FittedEdgs,FittedImage] = fitEdgs_f(image,edgPm,BObj,M,N)
%�������ƥ�䷽����������������ƽ������

%���У�
% FittedEdgs --- ������������ϲ��ֵ���Ĥ
% FittedImage --- ���������ɫ���ͼ�����
% image --- һ����ʾͼ��ľ���
% edgPm --- �������ߵ������ṹ�壬�ڴ�ֻʹ�����е�ǰ��������������������Ӵ��ĸ��˵����꼰�����������������������У�
% BObj --- һ�����������������ƥ�䷽��
% [M,N] --- ԭͼ��ߴ�

%���磺
% [edgPm,w,edgs,M,N,Q] = getEdge_f('planet01.png');
% BObj = optimizeByGA_f(edgPm,M,N);
% image = imread('planet01.png');
% [FittedEdgs,FittedImage] = fitEdgs_f(image,edgPm,BObj,M,N);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FittedEdgs = zeros(M,N);
image = double(image); %ת��ͼ����������
for i=1:BObj(end) %����ÿ��ƥ��������
    
    %��ʼ��
    xy1 = [edgPm(BObj(2*i-1)).Curve(:,2),edgPm(BObj(2*i-1)).Curve(:,1)]; %�ɵڶ˵�1�����߶ε������к������зֱ�õ�����1��x�������к�y��������
    xy2 = [edgPm(BObj(2*i)).Curve(:,2),edgPm(BObj(2*i)).Curve(:,1)]; %�ɵڶ˵�2�����߶ε������к������зֱ�õ�����2��x�������к�y��������
    p1 = [edgPm(BObj(2*i-1)).Beginning(2),edgPm(BObj(2*i-1)).Beginning(1)]; %�õ�һ�����𴦶˵㣨�˵�1���ĺ�������
    p2 = [edgPm(BObj(2*i)).Beginning(2),edgPm(BObj(2*i)).Beginning(1)]; %�õ���һ�����𴦶˵㣨�˵�2���ĺ�������
    
    %%������ɢ������
    figure(i)
    subplot(2,2,1)
    plot(xy1(:,1),xy1(:,2),'.',xy2(:,1),xy2(:,2),'.')
    title('Original discrete curve','FontSize',10,'FontAngle','italic');
    axis equal %ʹͼ��ĺ����������ʾ����һ��
    
    %�Ƕȹ�һ��
    v = [p2(1)-p1(1),p2(2)-p1(2)]; %���˵������
    angle = atan(v(2)/v(1)); %���˵�����������
    newXY1 = xy1*[cos(angle),-sin(angle);sin(angle),cos(angle)]; %�����������ת����1����������
    newXY2 = xy2*[cos(angle),-sin(angle);sin(angle),cos(angle)]; %�����������ת����2����������
    newp1 = newXY1(1,:); %�������ת�Ķ˵�1����
    newp2 = newXY2(1,:); %�������ת�Ķ˵�2����
    
    %��������1���ظ�ӳ�䲿��
    mean5XY1 = (newXY1(5:end,1)+newXY1(4:end-1,1)+newXY1(3:end-2,1)+newXY1(2:end-3,1)+newXY1(1:end-4,1))/5; %�õ���һ������x�����5�����
    dM5XY1 = mean5XY1(1:end-1) - mean5XY1(2:end); %����֡��õ�ÿ��������ֵ����һ��������ֵ�Ĳ�ֵ
    %position1 = find(dM5XY1(1)*dM5XY1<0); %�����������ֵ���Ÿı��λ��
    position1 = find(fix(10*dM5XY1)==0); %�����������ֵ��С��λ��
    if length(position1)~=0
        newXY1 = newXY1(1:position1(1)+4,:); %ֻ�����Ӷ˵㵽��һ��������5���ֵ�仯��С�ĵ�
    end
    
    %��������2���ظ�ӳ�䲿��
    mean5XY2 = (newXY2(5:end,1)+newXY2(4:end-1,1)+newXY2(3:end-2,1)+newXY2(2:end-3,1)+newXY2(1:end-4,1))/5; %�õ���һ������x�����5�����
    dM5XY2 = mean5XY2(1:end-1) - mean5XY2(2:end); %����֡��õ�ÿ��������ֵ����һ��������ֵ�Ĳ�ֵ
    %position2 = find(dM5XY2(1)*dM5XY2<0); %�����������ֵ���Ÿı��λ��
    position2 = find(fix(10*dM5XY2)==0); %�����������ֵ��С��λ��
    if length(position2)~=0
        newXY2 = newXY2(1:position2(1)+4,:); %ֻ�����Ӷ˵㵽��һ��������5���ֵ�仯��С�ĵ�
    end
     
    %%���ƽǶȹ�һ����ȥ��Ӱ�����ɢ������
    figure(i)
    subplot(2,2,2)
    plot(newXY1(:,1),newXY1(:,2),'.',newXY2(:,1),newXY2(:,2),'.')
    title('Unitary discrete curve','FontSize',10,'FontAngle','italic');
    axis equal %ʹͼ��ĺ����������ʾ����һ��
    
    %���
    f = fit([newXY1(:,1);newXY2(:,1)],[newXY1(:,2);newXY2(:,2)],'poly5'); 
    X = newp1(1):(newp2(1)-newp1(1))/(10*abs(newp2(1)-newp1(1))):newp2(1); %��x����Ҳȡ��С�����Ӷ�ʹ��ϳ���������һ���������ӷ���
    Y = f(X); %�����Ӧ��y����

    %%������ϲ��ֵ�����
    figure(i)
    subplot(2,2,3)
    plot(newXY1(:,1),newXY1(:,2),'.',newXY2(:,1),newXY2(:,2),'.')
    hold on
    plot(X,Y,'r.')
    title('Fitted curve','FontSize',10,'FontAngle','italic');
    axis equal %ʹͼ��ĺ����������ʾ����һ��
     
    %�ǶȻ�ԭ
    oXY = [X',Y]*[cos(angle),sin(angle);-sin(angle),cos(angle)]; %�õ��ǶȻ�ԭ������������������
    oXY = round(oXY); %ȡ��
     
    %%���ƽǶȻ�ԭ����������
    figure(i)
    subplot(2,2,4)
    plot(xy1(:,1),xy1(:,2),'.',xy2(:,1),xy2(:,2),'.')
    hold on
    plot(oXY(:,1),oXY(:,2),'r.')
    title('Recovered fitted curve','FontSize',10,'FontAngle','italic');
    axis equal %ʹͼ��ĺ����������ʾ����һ��
    
     
    %������������ͼ���е���Ĥ
    FittedEdgs(sub2ind(size(FittedEdgs),oXY(:,2),oXY(:,1))) = 1; %�������������λ�ñ�1
    
    %��ԭͼ��Ϊ�������ɫ
    reducedXY = [oXY(find((sum((abs(oXY(1:end-1,:)-oXY(2:end,:)))'))'~=0),:);oXY(end,:)]; %%
    %ȥ��ԭ����������������е��ظ����ꡣ...
    %�Ȳ�֣��ҵ���������һ���겻ͬ��������ţ�ԭ�����У���Щ���λ���ϵ����궼�ǻ��಻�ظ��ģ�������ԭ�����е����һ�����ꡣ
    for l = 1:size(image,3)
        tstep = (image(reducedXY(1,2),reducedXY(1,1),l) - image(reducedXY(end,2),reducedXY(end,1),l))...
            /(size(reducedXY,1)-1); %�õ��˴�ɫֵ�仯����ʱ����
        for j = 2:size(reducedXY,1)-1
            image(reducedXY(j,2),reducedXY(j,1),l) = image(reducedXY(j-1,2),reducedXY(j-1,1),l) - tstep; %������������ɫ
        end
    end
    
      
    %���������ԭ�����Ӵ�����Ԥ���Բ�ֵ
    %FittedEdgs(p1(2):(Y(1)-p1(2)+1)/abs(Y(1)-p1(2)+1):Y(1),X(1)) = 1;
    %FittedEdgs(p2(2):(Y(end)-p2(2)+1)/abs(Y(end)-p2(2)+1):Y(end),X(end)) = 1;
    %[�˴���������������ټ�1��Ϊ�˷�ֹ�����غ�ʱ���ֳ���0�����]
end
FittedEdgs = FittedEdgs==1; %2ֵ��
FittedImage = uint8(image); %ת����ͼ��������������

%%�������յ��������Ĥͼ
%imview(FittedEdgs)

%%���ƺ�����ߵ���Ƭͼ�� imview(uint8(image))
imshow(FittedImage)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     