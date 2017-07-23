function MeasPairs=SensorTask_GreedyTime_exhaust_jointcov(MeasPairs,Satellites,Radars,Constants,Tk,Tk1,method)
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

%%

% FIMG=zeros(size(XsigSat,1),Radmodel.Nrad,length(curk:fink));
MI=zeros(Constants.Nsat,Constants.Nrad);
MP=MeasPairs;


for k=Tk:Tk1
    disp(strcat( 'tasking for time step :',num2str(k) ) )
    MI=zeros(Constants.Nsat,Constants.Nrad);
    
    Satellites=Propagate_sattask(Satellites,Constants,k-1,k,'ut');
    
    for Ns=1:1:Constants.Nsat
        for Nr=1:Constants.Nrad
            mk=Satellites{Ns}.mu(k,:)';
            Pk=reshape( Satellites{Ns}.P(k,:),Satellites{Ns}.fn,Satellites{Ns}.fn );
            [x,w]=qd_pts(mk,Pk);
            H=zeros(size(x,1),1);
            for msi=1:1:size(x,1)
                [gg,hh]=Radars{Nr}.G( x(msi,:)', Radars{Nr}.PolarPositions, Radars{Nr}.hn, Radars{Nr}.ConeAngle,Radars{Nr}.MaxRange,Radars{Nr}.penalty);
                H(msi)=hh;
            end
            if sum(isnan(H))<length(H)/2
                MP{k}=zeros(Constants.Nsat,Constants.Nrad);
                MP{k}(Ns,Nr)=1;
                ymeas=NaN; % for pseudo update

                SS=MeasUpdate_sattask(MP,Satellites,Radars,Constants,k,ymeas,'ut','pseudoupdate');
                Pu=reshape( SS{Ns}.P(k,:),SS{Ns}.fn,SS{Ns}.fn );
                MI(Ns,Nr)=max(-0.5*log(det(Pu)/det(Pk)) ,0 ) ;
            else
                MI(Ns,Nr)=-1000000000;
            end
        end
    end
    % Now solve the binary integer program for this time step
    f=reshape(MI,1,Constants.Nsat * Constants.Nrad);
    ONZR=zeros(Constants.Nsat * Constants.Nrad,1);
    for i=1:Constants.Nsat
        ONZR(i)=1;
    end
    A=ONZR';
    for i=1:Constants.Nrad-1
        ONZR=circshift(ONZR,Constants.Nsat);
        A=vertcat(A,ONZR');
    end
    
    b=ones(Constants.Nrad,1);
    try
        mmu = bintprog(-f,A,b);
    catch
        keyboard
    end
    
    MeasPairs{k}=reshape(mmu,Constants.Nsat,Constants.Nrad);
    
    ymeas=NaN;
    Satellites=MeasUpdate_sattask(MeasPairs,Satellites,Radars,Constants,k,ymeas,'ut','pseudoupdate');
    
end



end