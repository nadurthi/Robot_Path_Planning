r=[7000,0,0,0,-1.03741,7.4771];
P=diag([0.01,0.01,0.01,1e-6,1e-6,1e-6]);
X=mvnrnd(r,P,1e5);
A=zeros(1e5,6);
for i=1:1:1e5
% A(i,:) = RVtoCOEs(X(i,1:3),X(i,4:6));
% [a,e,hh,omg,Omg,M]=cart2orb(X(i,1),X(i,2),X(i,3),X(i,4),X(i,5),X(i,6));
[a,e,hh,omg,Omg,M]= elorb(X(i,1:3),X(i,4:6));
A(i,:) =[a,e,hh,omg,Omg,M];
end
mean(A,1)
diag(sqrtm(cov(A)))