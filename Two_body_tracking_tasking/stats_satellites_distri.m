global tumin mu radiusearthkm xke j2 j3 j4 j3oj2  
   global opsmode
   opsmode= 'a';
global idebug dbgfile

endtime.year=2009;
endtime.month=1;
endtime.day=1;
endtime.hr=1;
endtime.min=01;
endtime.sec=01;
dtmin=0.5;

fid = fopen('cosmos_2251.txt');
X=zeros(2,6);
P=zeros(2,5);
k=1;
for i=1:1:9397
tline1 = fgets(fid);
if strcmp('1',tline1(1))
    tline2 = fgets(fid);
    tline2=strcat(tline2,'     0.00      4320.0        360.00');
    [satrec, startmfe, stopmfe, dtmin] = twoline2rv_modified(72,tline1,tline2,'m','e',endtime,dtmin);
    if  satrec.EpochDate(1)>=2008 && satrec.EpochDate(2)>=12 && satrec.EpochDate(1)<2009 
        [satrec, rr, vv] = sgp4(satrec,stopmfe);
        P(k,:)=[satrec.inclo,satrec.nodeo,satrec.ecco,satrec.a,satrec.argpo];
        X(k,:)=[rr(:)',vv(:)'];
        k=k+1;
    end
end

end
    
   diag(cov(P)) 
    
    fclose(fid);
    
%  Pir(i,Omg,e,a,omg)=[4.11600055847605e-011,0.00419439138141434,...
%       3.59640552995392e-012,8.55694161919788e-013,0.000769465668652812];   
% Pcm(i,Omg,e,a,omg)=[2.4762775560723e-010,0.0742804074381186,...
 %       1.43938233117484e-010,4.88521425304502e-014,0.088464833029821];   
% PrvCM=[1.50353408450567,2.98590210878029,17.4217801337083, ...
%     3.13134144653623e-007,2.00834441509321e-005,2.85465919597094e-006]  
% PrvIR=[65.4230995641526,43.9573041898352,93.8282101652225, ...
%     5.82803932257116e-005,5.6081579260967e-005,0.000109198004026865]  