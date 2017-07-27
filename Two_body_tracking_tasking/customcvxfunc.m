function y=customcvxfunc(mmu)

inv(diag(mmu))


real(mmu(1))
A=[0,1;-1,0];
A=eye(2)+0.5*A;
P=eye(2);
mux=[1;1];
iR=inv(diag([2^2,1^2]));
R=diag([2^2,1^2]);
h=@(x)[norm(x);x(1)];
hn=2;
H=eye(2);
mmu;
for i=1:1:length(mmu)
%     trace(inv(P))
    [X,w]=UT_sigmapoints(mux,P,2);
    [N,n]=size(X);
    for j=1:1:N
      X(j,:)=A*X(j,:)';
    end
    
      W=repmat(w,1,n);
      mux=sum(W.*X,1)';
      MU=repmat(mux',N,1);
      Xc=X-MU;
      P=Xc'*(W.*Xc);
      
      [X,w]=UT_sigmapoints(mux,P,2);
      Z=zeros(N,hn);
      for j=1:1:N
         Z(j,:)=h(X(j,:));
      end
      [Nz,nz]=size(Z);
      W=repmat(w,1,nz);
      muz=sum(W.*Z,1)';
      MUz=repmat(muz',Nz,1);
      Zc=Z-MUz;
      Pz=Zc'*(W.*Zc)+R;
      
      Pcc=0;
      for j=1:1:N
         Pcc=Pcc+w(j)*(X(j,:)'-mux(:))*(Z(j,:)-muz(:)'); 
      end
      K=Pcc/Pz;
      P=P-mmu(i)*K*Pz*K';
  
end
y=trace(inv(P));

end