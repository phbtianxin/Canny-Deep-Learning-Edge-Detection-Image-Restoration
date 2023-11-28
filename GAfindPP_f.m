function [ppc,Rpatch] = GAfindPP_f(image,w,edgs,FittedEdgs,mpc,Upatch,Dpatch,tp)
%���ҳ���Ӧ�ڸ���������ƥ������������
%ʹ���ⲿ������parallelMove()

%ǰ�����
%1����������ͼ��ı߽��㹻Զ
%2����ɨ��������ֻ�ٴ���������ߴ�һ�µ�����ƥ���

%���У�
% ppc --- һ��m��2�о������i�м�¼���Ƕ�Ӧ�ڵ�i����������ƥ������������
% Rpatch --- һ��m�����ÿһ����һ������ƥ����п���䲿�ֵ���Ĥ
% image --- һ����ʾͼ��ľ���
% w --- ��������Ĥ
% edgs --- ԭ�����������Ĥ
% FittedEdgs --- ����ϵ���������Ĥ
% mpc --- һ��������ʾ��ǰ�������������ȼ������������������е�m��2�о���
% Upatch --- һ��m�����ÿһ����һ�����������п��ò��ֵ���Ĥ
% Dpatch --- һ��m�����ÿһ����һ����������������䲿�ֵ���Ĥ
% tp --- �������ĵ�������һ���ߵľ���

%ʾ��:
% [edgPm,w,edgs,M,N,Q] = getEdge_f('planet01.png');
% BObj = optimizeByGA_f(edgPm,M,N);
% FittedEdgs = fitEdgs_f(edgPm,BObj,M,N);
% [mpc,Upatch,Dpatch,tp] = GAfindMPC_f(w,edgs,FittedEdgs,5);
% image = imread('planet01.png');
% [ppc,Rpatch] = GAfindPP_f(image,w,edgs,FittedEdgs,mpc,Upatch,Dpatch,tp);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tr = 8; %�������ĵ�ɨ������һ��������Ӧ��ɨ�����ĵľ��루����Ҫ��tp��

[M,N,L] = size(image);
uw = w|edgs|FittedEdgs; %����һ����������������������ߵ���Ĥ����������Ĥ��

ppc = []; %��ʼ�����ƥ����������������
for k = 1 : size(mpc,1) %���ڸ����ĸ�����
    
    %��ö��ڸ�������ĵ�ɨ�������ڵĿ��ò��ֵ���Ĥ
    %ǰ�᣺����������㹻�ָǰɨ������
    region = zeros(M,N); %��ʼ������ɨ������Ĥ��������������������������ߡ�����ͬһ���ӷ�����
    region(mpc(k,1)-tr-tp:mpc(k,1)+tr+tp,mpc(k,2)-tr-tp:mpc(k,2)+tr+tp) = 1; %���������ɨ������Ĥ
    region = region==1; %2ֵ��
    region = region&(~uw); %�ڿ�������Ĥ���ڳ���ɨ�����ڵĿ��ò���
    region(mpc(k,1),mpc(k,2)) = 1; %���ϸ���������
    region = bwlabel(region,4); %�����������򻮷�Ϊ�����ӷ���
    region = region==region(mpc(k,1),mpc(k,2)); %����������������������ӷ�����1
    region(mpc(k,1),mpc(k,2)) = 0; %Ĩȥ����������
    
    UDp = Upatch(:,:,k)|Dpatch(:,:,k); %�õ�������Ŀ��ò��ֺ�����䲿�ֵ�������Ĥ
    UDp = UDp(mpc(k,1)-tp:mpc(k,1)+tp,mpc(k,2)-tp:mpc(k,2)+tp); %�ٳ�������Ĥ�е����鲿��
    region = imerode(region,UDp); %�ø�����Ŀ��ò��ֺ�����䲿��һ��ʴ��ʼɨ�������õ���ȫ�ʺϸ������ɨ������
    regionSub = []; %���³�ʼ��
    [regionSub(:,1),regionSub(:,2)] = find(region==1); %��ȡ����ɨ������������

    minSV = inf; %%��ʼ����Сƽ�����ƽ����
    ppc(k,:) = [0,0]; %��ʼ����Ӧ�ڵ�i����������ƥ������������
    Rpatch(:,:,k) = zeros(M,N); %��ʼ����Ӧ�ڸ������ƥ����п���䲿�ֵ���Ĥ
    for l = 1:size(regionSub,1)
        i = regionSub(l,1); %������
        j = regionSub(l,2); %������
        
        %��������Ŀ��ò��ֺ�����䲿��ͬʱƽ�Ƶ���ǰɨ��λ��
        Upw = parallelMove(Upatch(:,:,k),[i-mpc(k,1),j-mpc(k,2)]); %�õ�ƽ�Ƹ������п��ò��ֵ���ǰɨ��λ�õ���Ĥ
        Dpw = parallelMove(Dpatch(:,:,k),[i-mpc(k,1),j-mpc(k,2)]); %�õ�ƽ�Ƹ�����������䲿�ֵ���ǰɨ��λ�õ���Ĥ
        Upw = Upw==1; %2ֵ��
        Dpw = Dpw==1; %2ֵ��
        
        if L==3
            Upw = cat(3,Upw,cat(3,Upw,Upw)); %��չ��������ͼ��ͬά
            oUpw = cat(3,Upatch(:,:,k),cat(3,Upatch(:,:,k),Upatch(:,:,k))); %��չ��������ͼ��ͬά
        else
            oUpw = Upatch(:,:,k); %�Ҷ�ͼ����
        end
        
        Up = Upw.*double(image); %������ͼ�����ڳ���ǰɨ��λ���еĿ��ò���
        Up = parallelMove(Up,[mpc(k,1)-i,mpc(k,2)-j]); %�������ڳ��Ĳ�����ƽ�ƻظ�����λ��
        oUp = oUpw.*double(image); %������ͼ�����ڳ����������еĿ��ò���
        SV = (sum(sum(sum((Up-oUp).^2))))/sum(sum(Upatch(:,:,k))); %���㵱ǰɨ�����������ͬ�����ò���֮���ƽ�����ƽ����
        
        if SV<minSV %�����ƽ�����ƽ����С�ڵ�ǰ��Сƽ�����ƽ����
            minSV = SV; %%���µ�ǰ��Сƽ�����ƽ����
            Rpatch(:,:,k) = Dpw; %���¶�Ӧ�ڸ������ƥ����п���䲿�ֵ���Ĥ
            ppc(k,:) = [i,j]; %���¶�Ӧ�ڵ�i����������ƥ������������
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
    