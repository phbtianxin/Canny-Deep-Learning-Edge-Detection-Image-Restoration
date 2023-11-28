function [mpc,Upatch,Dpatch,tp] = GAfindMPC_f(w,edgs,FittedEdgs)
%���ҳ���ǰ�������������ȼ�������������������
%��Ҫ��֤�ҳ�����Щ���ľ���������ȼ���
%��Ҫ��֤������������֮��ĺ�������������������2*tp��
%��Ҫ��֤ÿ���ҳ������ľ����ܵĶ�

%ǰ�����
%1����������ͼ��ı߽��㹻Զ

%���У�
% mpc --- һ��������ʾ��ǰ�������������ȼ������������������е�m��2�о���
% Upatch --- һ��m�����ÿһ����һ�����������п��ò��ֵ���Ĥ
% Dpatch --- һ��m�����ÿһ����һ����������������䲿�ֵ���Ĥ
% tp --- �������ĵ�һ���ߵľ���
% w --- ��������Ĥ
% edgs --- ԭ�����������Ĥ
% FittedEdgs --- ����ϵ���������Ĥ

%ʾ����
% [edgPm,w,edgs,M,N,Q] = getEdge_f('planet01.png');
% BObj = optimizeByGA_f(edgPm,M,N);
% FittedEdgs = fitEdgs_f(edgPm,BObj,M,N);
% [mpc,Upatch,Dpatch,tp] = GAfindMPC_f(w,edgs,FittedEdgs);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[M,N] = size(w);

%tp = 8; %�������ĵ�һ���ߵľ���
tp=4;

candidate = bwperim(w,8); %ȡ�������������������Ϊ��ѡ����
candidate = candidate&(~FittedEdgs); %����������ϵ����ز��ܱ���Ϊ��ѡ����
[canSub(:,1),canSub(:,2)]=find(candidate==1); %��ø���ѡ���ĵ�����

%%���ƺ�ѡ����ͼ
%imview(candidate)

mpc = []; %��ʼ���������������ȼ�������������������
mPRI = -1; %��ʼ�����������ȼ���ֵ����Ϊ-1��Ϊȷ��Upatchһ���ܱ���ʼ����
mpcW = zeros(M,N); %��ʼ��һ��ѡ���������ĵķֲ�ͼ��һ������һ�����㣩
for i=1:size(canSub,1)
    patch = ones(M,N); %��ʼ��һ��������Ĥ
    patch(canSub(i,1)-tp:canSub(i,1)+tp,canSub(i,2)-tp:canSub(i,2)+tp) = 0; %����i����ѡ�����ϵ�������0
    patch = patch==1; %2ֵ��
    patch = patch|w|edgs|FittedEdgs; %�ڳ��������г���������ĸ��������
    patch(canSub(i,1),canSub(i,2)) = 0; %������������
    patch = ~patch; %���ף�ʹ����������ֵΪ1
    patch = bwlabel(patch,4); %����������򻮷�Ϊ�����ӷ���
    patch = patch==patch(canSub(i,1),canSub(i,2)); %���������������������ӷ�����1
    patch(canSub(i,1),canSub(i,2)) = 0; %��ȥ��������
    PRI = sum(sum(patch)); %����������ȼ���������Ĥ�е�1ֵ���صĸ�������
    
    if PRI>mPRI %�������������ȼ��ȵ�ǰ���������ȼ���
        mPRI = PRI; %����µ�ǰ���������ȼ�
        mpc = canSub(i,:); %����mpc
        mpcW = zeros(M,N); %����ѡ���������ĵķֲ�ͼ
        mpcW(sub2ind(size(mpcW),mpc(:,1),mpc(:,2))) = 1; %����ѡ���������ĵķֲ�ͼ
        Upatch = []; %����Upatch
        Upatch(:,:,1) = patch; %����Upatch
    elseif (PRI==mPRI)&(sum(sum(mpcW(canSub(i,1)-2*tp:canSub(i,1)+2*tp,canSub(i,2)-2*tp:canSub(i,2)+2*tp)))==0) %%
        %����������е�ǰ������ȼ��Ҳ��������ִ��������ص�����
        mpc = [mpc;canSub(i,:)]; %������������������ӵ�����������������������
        mpcW(canSub(i,1),canSub(i,2)) = 1; %���ִ�ѡ���������ķֲ�ͼ����Ӹ����ĵ�
        Upatch = cat(3,Upatch,patch); %���������п��ò��ֵ���Ĥ������ӵ�Upatch��
    end
end
Upatch = Upatch==1; %2ֵ��

%��ø�����������䲿�ֵ���Ĥ
for i=1:size(mpc,1)
    tDpatch = zeros(M,N); %��ʼ��
    tDpatch(mpc(i,1)-tp:mpc(i,1)+tp,mpc(i,2)-tp:mpc(i,2)+tp) = 1; %����������������Ĥ
    tDpatch = tDpatch==1; %2ֵ��
    tDpatch = tDpatch&(~FittedEdgs); %����ϳ��������߽�������ָ�ɸ����ӷ�����������߾��������飩
    tDpatch = tDpatch&w; %�ڳ��������е�����䲿��
    tDpatch = bwlabel(tDpatch,4); %����Ϊ�����ӷ���
    tDpatch = tDpatch==tDpatch(mpc(i,1),mpc(i,2)); %ѡ�������������ڵ����ӷ���
    Dpatch(:,:,i) = tDpatch; %��¼
end

%%���Ƹ������п��ò���
%for i=1:size(Upatch,3)
%    imview(Upatch(:,:,i))
%    imview(Dpatch(:,:,i))
%end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    