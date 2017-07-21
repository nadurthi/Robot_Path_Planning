function err=RMSEMat(M1,M2)
a=M1-M2;
nr=size(a,1);
nc=size(a,2);

err=sqrt(sum(a.^2,2)/nc);
err=sqrt(sum(err.^2)/nr);


end