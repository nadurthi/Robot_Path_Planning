function [CovMaxTrace] = GetSatMetric(Satellites,Constants,Tk)



CovMaxTrace=zeros(1,length(1:Tk) );


for k=1:Tk
    S=ones(1,Constants.Nsat);
    for i=1:1:Constants.Nsat
        P=sqrtm( reshape( Satellites{i}.P(k,:) , Satellites{i}.fn, Satellites{i}.fn ) );
        S(i)=trace(P);
    end
    CovMaxTrace(k)=max(S);
end











end