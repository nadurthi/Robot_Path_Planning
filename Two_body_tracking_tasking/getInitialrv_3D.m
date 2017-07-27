function Xsat0=getInitialrv_3D(n,aL,aU)
i=1;
Xsat0=zeros(n,6);
while(i<=n)
    a = aL + (aU-aL).*rand(1,1);
    ecc = 0 + (0.9-0).*rand(1,1);
    inc=0 + (pi/2-0).*rand(1,1);
    Omega=0 + (2*pi-0).*rand(1,1);
    w=0 + (2*pi-0).*rand(1,1);
    nu=0;
    [r,v] = randv(a,ecc,inc,Omega,w,nu);
    
    if a*(1-ecc)>6500
        Xsat0(i,:)=[r',v'];
        i=i+1;
    end
end

end