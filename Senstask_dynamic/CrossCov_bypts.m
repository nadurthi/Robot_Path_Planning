function Pcc=CrossCov_bypts(x,z,w)
[N,nx]=size(x);
W=repmat(w,1,nx);
mx=sum(W.*x,1)';

[N,nz]=size(z);
W=repmat(w,1,nz);
mz=sum(W.*z,1)';


Pcc=0;
for i=1:1:N
    Pcc=Pcc+w(i)*(x(i,:)'-mx)*(z(i,:)'-mz)';
end


end