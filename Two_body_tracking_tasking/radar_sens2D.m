function h=radar_sens2D(x,Srad,Radmodel)
x=x(1:2);
x=x(:);
[xenu,~,~]=vec_radar_coordchange2D(x,Srad,Radmodel.RadPos,'ecef2local');
% x=x'-ecef_ref(Srad,:);
% xenu = ecef2enu(x*1e3,ecef_ref(Srad,:)*1e3)/1000;

r=norm(xenu);
th=atan2(xenu(1),xenu(2));


if abs(th)>Radmodel.SensParas(Srad,1) || r > Radmodel.SensParas(Srad,2)
    h=NaN;
    return;
end
if Radmodel.hn==2
    h=[r;th];
elseif Radmodel.hn==1
    h=r;
end

end