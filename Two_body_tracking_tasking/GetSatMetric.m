function [CovMaxTrace,RMSEpos,CovFrob] = GetSatMetric(Satellites,Constants,ytruth,Tk)



CovMaxTrace=zeros(length(1:Tk),3 );
RMSEpos=zeros(length(1:Tk),3 );
CovFrob=zeros(length(1:Tk),3 );

for k=1:Tk
    S0=ones(1,Constants.Nsat);
    S1=ones(1,Constants.Nsat);
    S2=ones(1,Constants.Nsat);
    
    for i=1:1:Constants.Nsat
        P=sqrtm( reshape( Satellites{i}.P(k,:) , Satellites{i}.fn, Satellites{i}.fn ) );
        S0(i)=trace(P);
        S1(i)=norm(P,'fro');
        
        S2(i)=norm(ytruth{i}(k,1:3)-Satellites{i}.mu(k,1:3) );
    end
    CovMaxTrace(k,:)=[min(S0),mean(S0),max(S0)];
    CovFrob(k,:)=[min(S1),mean(S1),max(S1)];
    
    RMSEpos(k,:)=[min(S2),mean(S2),max(S2)];
end











end