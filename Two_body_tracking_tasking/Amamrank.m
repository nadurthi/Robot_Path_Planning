Xprior=Data_TBP_at699cut8;
wprior=wcut8;
N=length(wprior);

y=[];
M=[];
ycell=cell(4,1);
Mcell=cell(4,1);
for i=1:1:4
[yy,MM]=Cal_moments_samples(Xprior,wprior,i,'central');
ycell{i}=yy;
Mcell{i}=MM;
y=vertcat(y,yy);
M=vertcat(M,MM);
end
y=[zeros(1,3);y];
M=[1;M];

A=zeros(size(y,1),N);
B=zeros(size(y,1),1);
% A(1,:)=ones(1,N);
% B(1)=1;

for i=1:1:size(y,1)
A(i,:)=prod(Xprior.^repmat(y(i,:),N,1),2);
B(i)=M(i);
end
