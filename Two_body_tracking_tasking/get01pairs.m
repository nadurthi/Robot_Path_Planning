function A=get01pairs(N)
A=zeros(2^N,N);
for i=0:1:2^N-1
    A(i+1,:)=de2bi(i,N);
    
end