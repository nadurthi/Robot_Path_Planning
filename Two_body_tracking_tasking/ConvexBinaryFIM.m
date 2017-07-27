% convex relaxation to solve FIM binary int

N=500;
FIM=2*randn(1,N).^2;
options = optimset('MaxFunEvals',30000);
mu=fmincon(@(x)-sum(x(:).*FIM(:))+0*sum(x(1:N/2)),[1,zeros(1,N-1)],[],[],[ones(1,floor(N/2)),zeros(1,floor(N/2))],1,zeros(1,N),ones(1,N),[],options);
% 1-exp(-1*x(:)))/(1-exp(-1))
mubin=bintprog(-FIM(:),[],[],[ones(1,floor(N/2)),zeros(1,floor(N/2))],1,[1,zeros(1,N-1)]);
plot(1:N,mu,'b+',1:N,mubin,'ro')

y=sin(linspace(0,2*pi,1000))+0.1*randn(1,1000);
D = diag(-1*ones(1000,1),0)+diag(1*ones(999,1),1);
K=2;
cvx_begin
    variable x(1000)
    minimize( sum_square(y'-x) )
    subject to
       norm(D*x,1) <= K;
cvx_end

x=0:0.01:1;
f=(1-exp(-5*x))/(1-exp(-5));

plot(x,f)