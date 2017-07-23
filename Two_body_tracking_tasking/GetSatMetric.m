function [CovMaxTrace] = GetSatMetric(Satellites,Constants)



CovMaxTrace=zeros(1,Constants.Ntimesteps);


for k=1:Constants.Ntimesteps
    S=ones(1,Constants.Nsat);
    for i=1:1:Constants.Nsat
        P=reshape( Satellites{i}.P(k,:) , Satellites{i}.fn, Satellites{i}.fn );
        S(i)=trace(P);
    end
    CovMaxTrace(k)=max(S);
end











end