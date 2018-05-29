% clc
close all
clear 
Nt=25;
Nmc=1000;
Ns=2;
dt=0.1;
Xk=zeros(Nt,Ns,Nmc);
Q=0.001*eye(2);
mu0=[1,1];
P0=0.01*eye(2);
Xmc0=mvnrnd(mu0,P0,Nmc);
Xk(1,:,:)=Xmc0';

for j=1:1:Nt-1
    for i=1:1:Nmc       
        Xk(j+1,:,i)=dummymodel(j,Xk(j,:,i),dt);
    end
    hnoise=zeros(1,Ns,Nmc);
    hnoise(1,:,:)=mvnrnd(zeros(1,2),Q,Nmc)';
    Xk(j+1,:,:)=Xk(j+1,:,:)+hnoise;
end
disp('done')
%%
figure
hold on
for i=1:Nmc/10
    plot(Xk(:,1,i),Xk(:,2,i),'ko--')
    plot(Xk(1,1,i),Xk(1,2,i),'ro')
    plot(Xk(end,1,i),Xk(end,2,i),'bo')
end
axis square
axis equal

%% Estimating the correlations of MC of time 5 and 15
t1=5;
t2=15;

X=zeros(Nmc,4);
for i=1:Nmc
    X(i,1:2)=Xk(t1,:,i);
    X(i,3:4)=Xk(t2,:,i);
end
mean(X)
cov(X)

%% Now estimating using quadratures




