function SMI = Compute_DynamicSens_MI_modf(currnr,abspath,Sensors,Targets,time,Tsteps,Tk)
% Compute the joint MI
% Now the code uses UT sigma points to estimate joint MI
% for everypath just regenerate all the joint sigma points from UT
global kappa
kappa=1;
qd_pts=@(m,P)UT_sigmapoints(m,P,2);

Ntseps=length(Tsteps);

for k=1:Ntseps
    Sensors.xf{currnr}(Tsteps(k),:)=abspath(k+1,:);
end
            


MI= zeros(1,Targets.Ntargs);

for Ntarg =1:Targets.Ntargs
    mk=Targets.xf{Ntarg}(Tk-1,:)';
    Pk=reshape(Targets.Pf{Ntarg}(Tk-1,:),Targets.fn(Ntarg),Targets.fn(Ntarg));
    
    % two stage: first compute for X
    mu_joint=[mk;zeros(Targets.fn(Ntarg)*(Ntseps),1)];
    
    statevars = zeros(Ntseps+1,2);
    statevars(1,:)=[1,Targets.fn(Ntarg)];
    for k=2:Ntseps+1
        statevars(k,:)=[statevars(k-1,2)+1,statevars(k-1,2)+Targets.fn(Ntarg)];
    end
    processnoises=zeros(Ntseps,2);
    processnoises(1,:)=[length(mk)+1,length(mk)+Targets.fn(Ntarg)];
    for k=2:Ntseps
        processnoises(k,:)=[processnoises(k-1,2)+1,processnoises(k-1,2)+Targets.fn(Ntarg)];
    end
    
    P_joint=Pk;
    for k=1:Ntseps
        P_joint=blkdiag(P_joint,Targets.Q{Ntarg});
    end
    % now add the sensors part
    measnoises=[];
    for nr= 1:Sensors.Nsens
        for k=1:Ntseps
            if Sensors.xf{nr}(Tsteps(k),1)>-1
                p=length(mu_joint);
                mu_joint=[mu_joint;zeros(Sensors.hn(nr),1)];
                P_joint=blkdiag(P_joint,Sensors.R{nr});
                measnoises=[measnoises;nr,k,p+1,p+Sensors.hn(nr)];
            end
        end
    end

    [D,W]=qd_pts(mu_joint,P_joint);

%     keyboard
    % D is the sigma points of the initial cond, process and meas noise
    % now simulate tge equations forward
    %     X=zeros(length(W),Targets.fn{Ntarg}*(Ntseps+1));
    X=zeros(size(D));
    G=zeros(Ntseps,Sensors.Nsens,length(W)); % G is used to check if target is visible to the sensor at time k
    for i=1:length(W)
        x=D(i,statevars(1,1):statevars(1,2))';
        X(i,statevars(1,1):statevars(1,2))= x';
        % propagate the state sequentially
        for k=2:Ntseps+1
            x=Targets.f{Ntarg}(time.dt,x)+D(i,processnoises(k-1,1):processnoises(k-1,2))';
            X(i,statevars(k,1):statevars(k,2))= x';
        end
        
        % Now get the measurement variables
        for s=1:size(measnoises,1)
            nr=measnoises(s,1);
            k=measnoises(s,2);
            a=measnoises(s,3);
            b=measnoises(s,4);
            x=X(i,statevars(k+1,1):statevars(k+1,2));
            [y,g,~] = Sensors.h{nr}( x, Sensors.xf{nr}(Tsteps(k),:), Sensors.FOV{nr} );
            if g==1
                G(k,nr,i)=G(k,nr,i)+1;
            end
            X(i,a:b)=y+D(i,a:b)';

        end
    end
    [mu,Pu]=MeanCov(X,W);
    % now remove the sensors that cannot see the target
    delcols=[];
    for s=1:size(measnoises,1)
        nr=measnoises(s,1);
        k=measnoises(s,2);
        a=measnoises(s,3);
        b=measnoises(s,4);
        c=sum(G(k,nr,:))/length(W);
        if c<0.25
            delcols=[delcols,a:b];
        else
            Pu(a:b,a:b)=Pu(a:b,a:b)+c*Sensors.R{nr}+(1-c)*Sensors.Rout{nr};
        end
    end
    
    Pu(:,delcols)=[];
    Pu(delcols,:)=[];
    pn=statevars(end,2);
    
    P=Pu(1:pn,1:pn);
    Z=Pu(pn+1:end,pn+1:end);
    S=Pu(1:pn,pn+1:end);
    if isempty(Z) || det(Z)==0 || det(Z- S'*(P\S) )==0
        MI(Ntarg)=0;
    else
        MI(Ntarg)=0.5*log( det(Z)/det(Z- S'*(P\S) ) );
    end
    
end

SMI = sum(MI);




