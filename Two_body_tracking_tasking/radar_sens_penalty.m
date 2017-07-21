function [g,h]=radar_sens_penalty(x,Srad,Radmodel,pen)
x=x(1:3);
x=x(:);
[xenu,~,~]=vec_radar_coordchange(x,Srad,Radmodel.RadPos,'ecef2local');
% x=x'-ecef_ref(Srad,:);
% xenu = ecef2enu(x*1e3,ecef_ref(Srad,:)*1e3)/1000;

r=norm(xenu);
th=atan2(xenu(1),xenu(2));
phi=atan2(sqrt(xenu(1)^2+xenu(2)^2),xenu(3));

if abs(phi)>Radmodel.SensParas(Srad,1) || r > Radmodel.SensParas(Srad,2)
    g=pen;
    h=NaN;
%     return;
else
    g=1;
    h=0;
end
% if Radmodel.hn==3
%     h=[r;th;phi];
% elseif Radmodel.hn==2
%     h=[r;th];
% elseif Radmodel.hn==1
%     h=r;
% end

end