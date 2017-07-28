function MeasPairs=SensorTask_GreedyTime_exhaust_jointcov(MeasPairs,Satellites,Radars,Constants,Tk,TkF,method)
%%
% - Tk-1 is fully updated time step.
% - Including Tk and Tk1
% - Over all time
% - Greedy over time
% - All independent at each time, solve binary integer problem
% - Next time step , use the pseudo updated covariances
% - Compute all mutual informations pairs for current time step conditioned
% on previous time step



% MeasPairs{1}(i,j)
%               j=1,rad1  j=2,rad2   j=3,rad3  j=4,rad4   ...
% i=1 sat1
% i=2 sat2
% i=3 sat3
% .
% .
% .


%%
opt = odeset('reltol',1e-12,'abstol',1e-12);
global kappa
kappa=1;
switch lower(method)
    case 'ut'
        qd_pts=@(m,P)UT_sigmapoints(m,P,2);
    case 'cut4'
        qd_pts=@conjugate_dir_gausspts;
    case 'cut6'
        qd_pts=@conjugate_dir_gausspts_till_6moment_scheme2;
    case 'cut8'
        qd_pts=@conjugate_dir_gausspts_till_8moment;
    case 'gh'
        qd_pts=@(m,P)GH_pts(m,P,para);
    otherwise
        error('smthg is wrong: DONT ask me what')
end

%% First get sall the covariances
% ZCov=cell(length(Tk:TkF),length(Tk:TkF),Constants.Nrad,Constants.Nrad,Constants.Nsat);
% TCov=cell(length(Tk:TkF),length(Tk:TkF),Constants.Nrad,Constants.Nsat);
% PCov=cell(length(Tk:TkF),length(Tk:TkF),Constants.Nsat);

Wsig=cell(1,Constants.Nsat);
MPnoviz=cell(1,Constants.Nsat);
Psig=cell(1,Constants.Nsat);
Zsig=cell(1,Constants.Nsat);
for i=1:Constants.Nsat
    Psig{i}=cell(1,length(Tk:TkF));
    Zsig{i}=cell(length(Tk:TkF),Constants.Nrad);
    MPnoviz{i}=cell(length(Tk:TkF),Constants.Nrad);
end
% Satellites_prop=Satellites;


MP=MeasPairs;
for k=Tk:TkF
   MeasPairs{k}=-1*ones(Constants.Nsat,Constants.Nrad); 
   MP{k}=-1*ones(Constants.Nsat,Constants.Nrad); 
end



%% First get all the satellite sigma points for all time steps
parfor i=1:1:Constants.Nsat
    i
    mk=Satellites{i}.mu(Tk,:)';
    Pk=reshape( Satellites{i}.P(Tk,:),Satellites{i}.fn,Satellites{i}.fn );
    
    F=Satellites{i}.f;
    for k=Tk:TkF
        if k==Tk
            [x,w]=qd_pts(mk,Pk);
            Psig{i}{1}=x;
            Wsig{i}=w;
        else
            Psig{i}{k-Tk+1}=zeros( length(Wsig{i}),Satellites{i}.fn );
            
            for msi=1:1:length(Wsig{i})
                [~,xx]=ode45(F ,Constants.Tvec([k-1,k]),Psig{i}{k-Tk}(msi,:)',opt);
                Psig{i}{k-Tk+1}(msi,:)=xx(end,:);
            end
        end
        
        % getting meas zig points
        
        for j=1:1:Constants.Nrad
            Zsig{i}{k-Tk+1,j}=NaN;
            ZZ=zeros( length(Wsig{i}), Radars{j}.hn );
            H=zeros( length(Wsig{i}), 1 );
            for msi=1:1:length(Wsig{i})
                ZZ(msi,:)=Radars{j}.h( Psig{i}{k-Tk+1}(msi,:)' , Radars{j}.PolarPositions, Radars{j}.hn );
                [gg,hh]=Radars{j}.G( Psig{i}{k-Tk+1}(msi,:)', Radars{j}.PolarPositions, Radars{j}.hn, Radars{j}.ConeAngle,Radars{j}.MaxRange,Radars{j}.penalty);
                H(msi)=hh;
                
            end
            Zsig{i}{k-Tk+1,j}=ZZ;
            if sum(isnan(H))>length(H)/2
                MPnoviz{i}{k,j}=0;
            end
        end
        
%         % cross cov of state with previous time state
%         for pk=Tk:k
%             if k==pk
%                 %(Xt1,Xt2,sat)
%                 PCov{pk-Tk+1,k-Tk+1,i}=CrossCov_bypts(Psig{pk-Tk+1,i}.X,Psig{k-Tk+1,i}.X,Psig{1,i}.W)+Satellites{i}.Q;
%             else
%                 PCov{pk-Tk+1,k-Tk+1,i}=CrossCov_bypts(Psig{pk-Tk+1,i}.X,Psig{k-Tk+1,i}.X,Psig{1,i}.W);
%             end
%             %             PCov{k-Tk+1,pk-Tk+1,i}=PCov{pk-Tk+1,k-Tk+1,i}';
%         end
%         
%         % cross cov of meas with previous time state
%         for jk=1:1:Constants.Nrad
%             for pk=Tk:k
%                 
%                 for jpk=1:1:Constants.Nrad
%                     if pk==k && jpk>jk
%                         continue
%                     end
%                     
%                     if (k==pk) && (jk==jpk)
%                         % (Ytk-1,Ytk,rad1,rad2,sat)
%                         ZCov{pk-Tk+1,k-Tk+1,jpk,jk,i}=CrossCov_bypts(Zsig{pk-Tk+1,jpk,i}.Z,Zsig{k-Tk+1,jk,i}.Z,Psig{1,i}.W)+Radars{jk}.R ;
%                     else
%                         ZCov{pk-Tk+1,k-Tk+1,jpk,jk,i}=CrossCov_bypts(Zsig{pk-Tk+1,jpk,i}.Z,Zsig{k-Tk+1,jk,i}.Z,Psig{1,i}.W);
%                     end
%                     %                     ZCov{k-Tk+1,pk-Tk+1,j1,j2,i}=ZCov{pk-Tk+1,k-Tk+1,j1,j2,i}';
%                 end
%             end
%         end
%         
%         % cross cov of state and meas all rad all prev time meas
%         for pk=Tk:k
%             for j=1:1:Constants.Nrad
%                 % (Xt,Yt,rad,sat)
%                 TCov{pk-Tk+1,k-Tk+1,j,i}=CrossCov_bypts(Psig{pk-Tk+1,i}.X,Zsig{k-Tk+1,j,i}.Z,Psig{1,i}.W);
%                 %                 TCov{pk-Tk+1,k-Tk+1,j,i}=
%             end
%         end
        
    end
    
    
end

for i=1:1:Constants.Nsat
    for k=Tk:TkF
        for j=1:1:Constants.Nrad 
            if isempty(MPnoviz{i}{k,j})==0
                MP{k}(i,j)=MPnoviz{i}{k,j};
            end
        end
    end
end

% keyboard

%% Greedy  in targets
% TraceCov=zeros(1,Constants.Nsat);
% for i=1:Constants.Nsat
%     Pk=reshape( Satellites{i}.P(Tk,:),Satellites{i}.fn,Satellites{i}.fn );
%     TraceCov(i)= trace(Pk);
% end
% [~,IDprioritySat]=sort(TraceCov,2,'descend');

II=1:Constants.Nsat*Constants.Nrad;
% Ps=[];
% Ts=[];
% Zs=[];
% for i=IDprioritySat
%     for pk=Tk:TkF
%         PP=[];
%
%         for k=Tk:TkF
%             if pk<=k
%                 PP=horzcat(PP,PCov{pk-Tk+1,k-Tk+1,i} );
%             else
%                 PP=horzcat(PP,PCov{k-Tk+1,pk-Tk+1,i}' );
%             end
%         end
%         Ps=vertcat(Ps,PP);
%
%         ZZ=[];
%         for jpk=1:Constants.Nrad
%             for k=Tk:TkF
%                 for jk=1:Constants.Nrad
%                     if pk<=k && jpk<=jk
%                         ZZ=horzcat(ZZ, ZCov{pk-Tk+1,k-Tk+1,jpk,jk,i} );
%                     else
%                         ZZ=horzcat(ZZ, ZCov{pk-Tk+1,k-Tk+1,jpk,jk,i}' );
%                     end
%                 end
%             end
%         end
%
%         Zs=vertcat(Zs,ZZ);
%     end
%
%
% end

%% Method by only taking sigma points
% Greedy by sensor and Greedy in time

III=eye(Constants.Nsat);
II=cell(1,Constants.Nsat);
for i=1:Constants.Nsat
   II{i}=III(:,i); 
end

for j=1:Constants.Nrad
%     MPrt=[];
%     for k=Tk:TkF
%         MPrt=vertcat(MPrt,MP{k}(:,j));
%     end
%     MPrt_vec=MPrt';
%     Nvar=MPrt_vec(MPrt_vec==-1);
    
    


    for k=Tk:TkF
        MMP=cell(1,Constants.Nsat);
        for i=1:Constants.Nsat
            MMP{i}=MeasPairs;
        end
        
        MI=zeros(1,Constants.Nsat);
        parfor i=1:Constants.Nsat
            MMP{i}{k}(:,j)=II{i};
            MMP{i}{k}( MP{k}(:,j)==0,j )=0;
            
            % compute the MI for MMP
%             keyboard
            MI(i)=ComputeJointMI(Satellites,Radars,Constants,MMP{i},Psig,Zsig,Tk,TkF,Wsig);

        end
        
        [Y,ind]=max(MI);
        
%         [j,k]
%         [Y,ind]
%         MI
%         disp('-++++++++++---')

        if Y==0 || isnan(Y)
            MeasPairs{k}(:,j)=0;
        else
            MeasPairs{k}(:,j)=0;
            if MP{k}(ind,j)==0
                MeasPairs{k}(ind,j)=0;
            else
                MeasPairs{k}(ind,j)=1;
            end
        end
        
        
    end  
    
end
% keyboard


end