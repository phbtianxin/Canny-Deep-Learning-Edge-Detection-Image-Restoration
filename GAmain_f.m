function RImage = GAmain_f(Image)
%��һ��ͼ�������������ɴ��ڱ�ʾ�����л����Ŵ��㷨���޸�

%���У�
% RImage --- һ����ʾͼ��ľ���
% Image --- һ���ַ���������һ��ͼ��

%ʾ����
% RImage = GAmain_f('planet01.png');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%image = imread(Image);

%%����ԭͼ��
%imview(image)
image=Image;
imshow(image);%function imview has been removed 

%��ȡ����������ͼ���������������нӴ��ĸ������������м���Ӧ�������
[edgPm,w,edgs,M,N,Q] = getEdge_f(Image);

%���ݲ�����Եĸ������ߵĲ���ͨ��GA�㷨�������ƥ�䷽��
BObj = optimizeByGA_f(edgPm,M,N);

%�������ƥ�䷽����������������ƽ������
[FittedEdgs,FittedImage] = fitEdgs_f(image,edgPm,BObj,M,N);

RImage = double(FittedImage); %��ʼ�����ͼ�����
%avoid endless loop
Record=[(sum(sum(w&~FittedEdgs))),1];
while sum(sum(w&~FittedEdgs))~=0 %ֻҪ����ϵ��������⻹����������
    
    %���ҳ���ǰ�������������ȼ�������������������
    [mpc,Upatch,Dpatch,tp] = GAfindMPC_f(w,edgs,FittedEdgs);
    
    %���ҳ���Ӧ�ڸ���������ƥ������������
    [ppc,Rpatch] = GAfindPP_f(RImage,w,edgs,FittedEdgs,mpc,Upatch,Dpatch,tp);
    
    for i=1:size(ppc,1)
        Rpw = Rpatch(:,:,i); %�õ�Ŀ��ƥ����еĿ���䲿�ֵ���Ĥ
        
        if Q==3 %����ǲ�ɫͼ��
            Rpw = cat(3,Rpw,cat(3,Rpw,Rpw)); %��չ��������ͼ��ͬά
        end
        
        Rp = Rpw.*RImage; %������ͼ�����ڳ����������������Ĳ���
        Rp = parallelMove(Rp,[mpc(i,1)-ppc(i,1),mpc(i,2)-ppc(i,2)]); %�������ڳ��Ĳ�����ƽ�ƻ�������λ��
        RImage = RImage + Rp; %���ͼ��
        w = w&(~(Rp(:,:,1)>0));%�����������Ĥ
    end
    %avoid endless loop
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
RImage = uint8(RImage); %��ԭ��ͼ������ݸ�ʽ

%%�����޸����ͼ�� imview(uint8(RImage))
%imview(RImage)
imshow(RImage);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        