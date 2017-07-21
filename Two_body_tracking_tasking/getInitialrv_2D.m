function Xsat0=getInitialrv_2D(n)
i=1;
Xsat0=zeros(n,4);
while(i<=n)
a = 14000 + (30000-14000).*rand(1,1);
ecc = 0 + (0.9-0).*rand(1,1);
inc=0;
Omega=0 + (2*pi-0).*rand(1,1);
w=0 + (2*pi-0).*rand(1,1);
nu=0;
    [r,v] = randv(a,ecc,inc,Omega,w,nu);

if a*(1-ecc)>6500
    Xsat0(i,:)=[r(1:2)',v(1:2)'];
i=i+1;
end
end

end