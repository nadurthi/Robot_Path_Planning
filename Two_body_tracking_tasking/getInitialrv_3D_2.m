function Xsat0=getInitialrv_3D_2(n)


% [a0,e0,i0,omg0,Omg0,M0]=cart2orb(7000,0,0,0,-1.0374090357,7.4771288355);
mue = 398601.2;

a0=8.1057e3;
e0=0.1364;
M0=0;
w0=0;
i0=1.7;
Omega0=0;

i=1;
Xsat0=zeros(n,6);
while(i<=n)
a = (a0-30) + ((a0+30)-(a0-30)).*rand(1,1);
ecc =  (e0-0.2) + ((e0+0.2)-(e0-0.2)).*rand(1,1);
inc= 0.5*i0 + (1.0*i0-0.5*i0).*rand(1,1);
Omega= 0.5*Omega0 + (1.0*Omega0-0.5*Omega0).*rand(1,1);
w=0.99*w0 + (1.01*w0-0.99*w0).*rand(1,1);
M=0 + (pi/20-0).*randn(1,1);
%     [r,v] = randv(a,ecc,inc,Omega,w,nu);
E=kepler(ecc,M);
% [r, v] = OE2XYZ(a, e0, E0sat, w0, inc(j), Om0, mue]);
XX= OE2XYZ([a, ecc, E, w, inc, Omega,mue]);
r=XX(1:3);
v=XX(4:6);

if a*(1-ecc)>6500
    Xsat0(i,:)=[r',v'];
i=i+1;
end
end

end