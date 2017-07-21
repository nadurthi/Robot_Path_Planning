function y=propagate_sigma_pts_2body(x,tT)

[N,n]=size(x);

if N<15
    y=zeros(size(x));
    for sdd=1:1:N
    [tt,xx]=ode45(@twoBody,tT,x(sdd,:)',opt);
%     [tt,xx]=twoBodyKeplerProp(tT,x(sdd,:)');
    y(sdd,:)=xx(end,:);
    end
else
    y=zeros(size(x));
     parfor sdd=1:N
    [tt,xx]=ode45(@twoBody,tT,x(sdd,:)',opt);
%     [tt,xx]=twoBodyKeplerProp(tT,x(sdd,:)');
    y(sdd,:)=xx(end,:);
    end   






end