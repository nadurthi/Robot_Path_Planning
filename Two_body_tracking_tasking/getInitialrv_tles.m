% %% From the following TLE data (latest before the collision)
% % 2/09/2009 at 11:57:36.8904960001419
% CML1='1 22675U 93036A   09040.49834364 -.00000001  00000-0  95251-5 0  7411';
% CML2='2 22675  74.0355  19.4646 0016027  98.7014 261.5952 14.31135643817415     0.00      4320.0        360.00';
% 
% % 2/09/2009 at 18:49:39.2819519997465  
% IRL1='1 24946U 97051C   09040.78448243  .00000153  00000-0  47668-4 0  4775';
% IRL2='2 24946  86.3994 121.7028 0002288  85.1644 274.9812 14.34219863597336     0.00      4320.0        360.00' ;
function X0=getInitialrv_tles()

global tumin mu radiusearthkm xke j2 j3 j4 j3oj2  
   global opsmode
   opsmode= 'a';
global idebug dbgfile

endtime.year=2009;
endtime.month=2;
endtime.day=10;
endtime.hr=16;
endtime.min=01;
endtime.sec=01;
dtmin=0.5;
fid = fopen('egtles.m', 'r');
tline=1;
X0=[];
while 1
tline = fgets(fid);
if tline==-1
    break
end
if strcmp(tline(1),'1')==1 && strcmp(tline(2),' ')==1
    L1=tline;
    L2=fgets(fid);
    L2=strcat(L2,'     0.00      4320.0        360.00');
[satrec, startmfe, stopmfe, dtmin] = twoline2rv_modified(72,L1,L2,'m','e',endtime,dtmin);
 [satrec, r, v] = sgp4(satrec,0);
 X0=vertcat(X0,[r(:)',v(:)']);
end

end
