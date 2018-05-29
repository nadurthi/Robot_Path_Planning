function MeasPairs=SensorTask_GreedyTime_exhaust_jointcov(MeasPairs,Satellites,Radars,Constants,Tk,TkF,method)
%%
% - get the tasking for Tk:TkF
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

Wsig=cell(1,Constants.Nsat);
% MPnoviz=cell(1,Constants.Nsat);
Psig=cell(1,Constants.Nsat);
Zsig=cell(1,Constants.Nsat);

MPnoviz=cell(1,Constants.Nsat);
for i=1:Constants.Nsat
    Psig{i}=cell(1,length(Tk:TkF));
    Zsig{i}=cell(length(Tk:TkF),Constants.Nrad);
    MPnoviz{i}=-1*ones(length(Tk:TkF),Constants.Nrad);
end
% Satellites_prop=Satellites;

Nttask=length(Tk:TkF);
Tvectask = Tk:TkF;

MP=MeasPairs;
for k=Tk:TkF
   MeasPairs{k}=-1*ones(Constants.Nsat,Constants.Nrad); 
   MP{k}=-1*ones(Constants.Nsat,Constants.Nrad); 
end

%% First get all the satellite sigma points for all time steps


for i=1:1:Constants.Nsat
    i
    mk=Satellites{i}.mu(Tk,:)';
    Pk=reshape( Satellites{i}.P(Tk,:),Satellites{i}.fn,Satellites{i}.fn )+Satellites{i}.Qmeas;
    
    F=Satellites{i}.f;
    for k=1:Nttask % Tk:TkF
        if k==1
            [x,w]=qd_pts(mk,Pk);
            Psig{i}{1}=x;
            Wsig{i}=w;
        else
            Psig{i}{k}=zeros( length(Wsig{i}),Satellites{i}.fn );
            
            for msi=1:1:length(Wsig{i})
                [~,xx]=ode45(F ,Constants.Tvec([Tvectask(k-1),Tvectask(k)]),Psig{i}{k-1}(msi,:)',opt);
                Psig{i}{k}(msi,:)=xx(end,:);
            end
        end
        
        % getting meas zig points
        try
        for j=1:1:Constants.Nrad
            Zsig{i}{k,j}=NaN;
            ZZ=zeros( length(Wsig{i}), Radars{j}.hn );
            H=zeros( length(Wsig{i}), 1 );
            for msi=1:1:length(Wsig{i})
                ZZ(msi,:)=Radars{j}.h( Psig{i}{k}(msi,:)' , Radars{j}.PolarPositions, Radars{j}.hn );
                [gg,hh]=Radars{j}.G( Psig{i}{k}(msi,:)', Radars{j}.PolarPositions, Radars{j}.hn, Radars{j}.ConeAngle,Radars{j}.MaxRange,Radars{j}.penalty);
                H(msi)=hh;  
            end
            Zsig{i}{k,j}=ZZ;
            if sum(isnan(H))>length(H)/2
                MPnoviz{i}(k,j)=0;
            end
        end
        catch
            keyboard
        end
        
    end
    
    
end

for k=1:Nttask %Tk:TkF
    for i=1:1:Constants.Nsat
        for j=1:1:Constants.Nrad 
            if MPnoviz{i}(k,j)==0
                MP{Tvectask(k)}(i,j)=MPnoviz{i}(k,j);
            end
        end
    end
end

% keyboard


%% Method by only taking sigma points
% Greedy by sensor

III=eye(Constants.Nsat);
II=cell(1,Constants.Nsat);
for i=1:Constants.Nsat
   II{i}=III(:,i); 
end

% Nttask=length(Tk:TkF);
% Tvectask

% greedy by sensor
% for each sensor do over all time stesp and all targets
for nrad=1:Constants.Nrad
%     for k=1:Nttask %Tk:TkF
        
        % generate all the possibilities
        
        % make a copy of the most update measParirs
%         MMP=MeasPairs;
        
        % generate all possibilities of size Constants.Nsat,Nttask
        % with a constraint that sum along rows is 1 for every column
        
        
%         MI=zeros(1,(Constants.Nsat)^Nttask);
        MImax=0;
        SenstaskMImax=[];
        A=zeros(1,Nttask);
        cnt=1;
        while(1)
            cnt=cnt+1;
            MMP=MeasPairs;
            
            % if flg ==1 then it is the end
            [A,flg] = Increment_index(A,Constants.Nsat);
            skipit=false;
            
            for k=1:Nttask
               if A(k)>0 && MP{Tvectask(k)}(A(k),nrad)==0
                   skipit=true;
                   break 
               end
            
               if A(k)==0 
                   MMP{Tvectask(k)}(:,nrad)=0;
               else
                   MMP{Tvectask(k)}(:,nrad)=0;
                   MMP{Tvectask(k)}(A(k),nrad)=1; % take care of no visibility
               end
               
            end
            
            if skipit==true
                if flg==1
                    break
                end
                continue
            end
            
            % MMP actually has the correct time steps
            % Psig and Zsig indexed from 1 for time
            MI=ComputeJointMI(Satellites,Radars,Constants,MMP,Psig,Zsig,Tk,TkF,Wsig);
            if MI > MImax
                MImax = MI;
                SenstaskMImax = MMP;
            end
            
            if flg==1
                break
            end
        end
        
        for k=1:Nttask
            MeasPairs{Tvectask(k)}(:,nrad)=0;
            if isempty(SenstaskMImax)==0
                MeasPairs{Tvectask(k)}(:,nrad) = SenstaskMImax{Tvectask(k)}(:,nrad);
            end
        end

        
        
%     end  
    
end
% keyboard


end