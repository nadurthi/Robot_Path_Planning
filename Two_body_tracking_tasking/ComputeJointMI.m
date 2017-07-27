function SumMI=ComputeJointMI(Satellites,Radars,Constants,MeasPairs,Psig,Zsig,Tk,TkF,Wsig)


MI=zeros(1,Constants.Nsat);


for i=1:1:Constants.Nsat
    X=[];
    pn=0;
    Q=[];
    R=[];
    
    XZ=[];
    for k=Tk:TkF
        for j=1:1:Constants.Nrad
            if MeasPairs{k}(i,j)==1
                XZ=horzcat(XZ,Zsig{i}{k-Tk+1,j});
                R=blkdiag(Q,Radars{j}.R);
            end
        end
        
    end
    if isempty(XZ)==1
        MI(i)=0;
        continue
    end
    
    for k=Tk:TkF
        [a,b]=size(Psig{i}{pk-Tk+1});
        X=horzcat(X,Psig{i}{pk-Tk+1});
        Q=blkdiag(Q,Satellites{i}.Q);
        pn=pn+b;
    end
    
    X=horzcat(X,XZ);
    
    
    
    if size(X,2)==pn
        MI(i)=0;
    else
        [mu,Pu]=MeanCov(X,Wsig{i});
        
        P=Pu(1:pn,1:pn)+Q;
        Z=Pu(pn+1:end,pn+1:end)+R;
        S=Pu(1:pn,pn+1:end);
        MI(i)=0.5*log( det(Z)/det(Z- S'*P*S ) );
    end
    
    
    
end


SumMI=sum(MI);


















end