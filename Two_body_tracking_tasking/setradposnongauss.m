Re=6378.1;

P0=blkdiag(0.01,0.01,0.01,1e-7,1e-7,1e-7);
Xsat0= [7200 0 0 1.0374090357 -1.3374090357 7.4771288355];
Xsat0=mvnrnd(Xsat0,P0,1000);
opt = odeset('reltol',1e-12,'abstol',1e-12);


Nsat=size(Xsat0,1);

ytruth=cell(Nsat,1);

tf=1*24*60*60; % 48 hours 
dt=1*60; % 10 mins time step
t0=0;
Tvec=t0:dt:tf;
nt=length(t0:dt:tf);

for i=1:Nsat
    
[tt,xx]=ode45(@twoBody,Tvec,Xsat0(i,:)',opt);
ytruth{i}=xx;
 i

end

 

