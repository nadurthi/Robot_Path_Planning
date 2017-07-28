function [v,R,p]=vec_radar_coordchange2D(v,nrad,RadPos,meth)
% Nrad=size(RadPos,1);
% for i=1:1:Nrad
% v=[v(:);0];
i=nrad;
th=RadPos(i,1);
r=RadPos(i,2);
x=r*cos(th);
y=r*sin(th);
p=[x;y];

R=[cos(th-pi/2),-sin(th-pi/2);sin(th-pi/2),cos(th-pi/2)]';
% zv=[x,y,0]/norm([x,y,0]);
% yv=-cross(zv,[0,0,1]);
% yv=yv/norm(yv);
% xv=-cross(zv,yv);
% xv=xv/norm(xv);
% 
% R=[xv;yv;zv];


% keyboard

if strcmp(meth,'ecef2local')
    v=R*(v(:)-[x;y]);
else
    v=R'*v(:)+[x;y];
    
end




% end
