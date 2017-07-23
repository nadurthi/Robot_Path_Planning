function Xsat0=getInitialrv_3D_2(n)


% [a0,e0,i0,omg0,Omg0,M0]=cart2orb(7000,0,0,0,-1.0374090357,7.4771288355);
a0=8.1057e3;
e0=0.1364;
M0=0;
w0=0;
i0=1.7;
Omega0=0;

i=1;
Xsat0=zeros(n,6);
while(i<=n)
a = (a0-10) + ((a0+10)-(a0-10)).*rand(1,1);
ecc =  (e0-0.05) + ((e0+0.05)-(e0-0.05)).*rand(1,1);
inc= 0.99*i0 + (1.01*i0-0.99*i0).*rand(1,1);
Omega= 0.99*Omega0 + (1.01*Omega0-0.99*Omega0).*rand(1,1);
w=0.99*w0 + (1.01*w0-0.99*w0).*rand(1,1);
M=0 + (pi/5-0).*rand(1,1);
%     [r,v] = randv(a,ecc,inc,Omega,w,nu);
E=kepler(M, ecc, 0, 1e-5);
[r, v] = OE2XYZ(a, ecc, E, w, inc, Omega);

if a*(1-ecc)>6500
    Xsat0(i,:)=[r',v'];
i=i+1;
end
end

end