function phik_k1=STMcovPROP(muk,Tveckk1)
muconst=398601.2;

%page 129 130 battin 
muk=muk(:)';

r0=muk(1:3);
v0=muk(4:6);
nr0=norm(r0);
nv0=norm(v0);

h=cross(r0,v0);
nh=norm(h);
p=nh^2/muconst;

a=(2/nr0-nv0^2/muconst)^(-1);
e=sqrt(1-p/a);
sig0=dot(r0,v0)/sqrt(muconst);

% r=p*nr0/(nr0+(p-nr0)*cos(th)-sqrt(p)*sig0*sin(th));
% 
% F=1-r/p*(1-cos(th));
% G=r*nr0/sqrt(muconst*p)*sin(th);
% Ft=sqrt(muconst)/(nr0*p)*(sig0*(1-cos(th)-sqrt(p)*sin(th)));
% Gt=1-nr0/p*(1-cos(th));
options = odeset('RelTol',1e-8,'AbsTol',1e-8);

[t,phik_k1]=ode45(@(t,phi)STMdiff(t,phi,Tveckk1(1),muk,muconst),Tveckk1,reshape(eye(6),36,1));

% phik_k1=reshape(phik1(end,:),6,6);

end


function dphi=STMdiff(t,phi,t0,x0,muconst)

options = odeset('RelTol',1e-12,'AbsTol',1e-12);
if t-t0==0
    x2=x0;
else
[t2,x2]=ode45(@(t,x)[x(4);x(5);x(6);-muconst/norm(x(1:3))^3*[x(1);x(2);x(3)]],[t0,t],x0,options);
end

nr=norm(x2(end,1:3));
r=x2(end,1:3)';

G=muconst/nr^5*(-nr^2*eye(3)+3*(r*r'));
dphi=[zeros(3,3),eye(3);G,zeros(3,3)]*reshape(phi,6,6);

dphi=reshape(dphi,36,1);


end