function C=Cofactor(P,n,ni,nj)

M=det(P([1:(ni-1),(ni+1):n],[1:(nj-1),(nj+1):n]));
C=(-1)^(ni+nj)*M;
end