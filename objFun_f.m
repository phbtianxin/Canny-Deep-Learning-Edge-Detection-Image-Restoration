function ObjV = objFun_f(Chrom,edgPm)

% N = 5;
% Chrom = crtbp(7, [N N-1 N-2 N-3 floor(N/2)]);
% edgPm = rand(5,4);
% ObjV = objFun_f(Chrom,edgPm);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T = 0.3; 
[Nind,Lind] = size(Chrom);

Chrom = Chrom+1;
Phenotype = 0*Chrom; 
for i=1:Nind
    tem = 1:Lind; 
    for j=1:2*Chrom(i,Lind) 
        Phenotype(i,j) = tem(Chrom(i,j)); 
        tem(Chrom(i,j)) = [];
    end
end
Nc = Chrom(:,Lind); 
        
ObjV = zeros(Nind,1);
k = 0;
while(sum(Nc)~=0)
    R = 0*Phenotype; 
    for i=1:Nind
        if Nc(i,1) ~= 0
            R(i,Phenotype(i,k+1))=1; 
            R(i,Phenotype(i,k+2))=-1; 
            Nc(i,1) = Nc(i,1)-1;
        end
    end
    ObjV = ObjV+(sum((abs(R*edgPm))'))'; 
    k = k+2;
end

ObjV = (T+ObjV)./(Chrom(:,Lind));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

