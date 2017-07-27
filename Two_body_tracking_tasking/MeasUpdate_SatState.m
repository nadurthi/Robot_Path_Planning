function SatState=MeasUpdate_SatState(Nsats,SatState,MeasPairs,Radmodel,Tk,Tvec,method,further)
% Nsats is the index of the satellites to be updated
% simply do the measurement update for all the satellites at the time step
% Tk

% make the measurement if ture, or just use the pseudo measurements
% Tk is absolute time in the Tvec

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


Nsats=length(SatState);

%%
if paras.TaskAll==true
    % overite teh tasking to task all sensors to all the satellites
    MeasPairs{Tk}=[];
    for Ns=Nsats
        MeasPairs{Tk}=vertcat(MeasPairs{Tk},[Ns*ones(Radmodel.Nrad,1) ,[1:Radmodel.Nrad]']);
    end
end

%%
for Ns=Nsats
    
    mk=SatState{Ns}.mu(Tk,:)';
    Pk=reshape(SatState{Ns}.P(Tk,:),6,6);
    
    
    [x,w]=qd_pts(mk,Pk);
    N=size(x,1);
    
    if sum(size(MeasPairs{Tk}))==0
        ind=[];
    else
        ind=find(MeasPairs{Tk}(:,1)==Ns);  %[satellites,sensors/radars]
    end
    
    
    
    if isempty(ind)==0
        Srad=MeasPairs{Tk}(ind,2);
        %         Y=zeros(size(x,1),Radmodel.hn*length(Srad));
        Y=[];
        ym=[];
        RR=[];
        GG=[];
        for nr=1:1:length(Srad)
            ZZ=zeros(N,Radmodel.hn);
            G=zeros(N,1);
            H=zeros(N,1);
            for msi=1:1:N
                if isreal(x(msi,:))==0
                    keyboard
                end
                ZZ(msi,:)=Radmodel.h(x(msi,:)',Srad(nr));        
                [gg,hh]=Radmodel.G(x(msi,:)',Srad(nr));
                G(msi)=gg;
                H(msi)=hh;

            end
            Y=horzcat(Y,ZZ);
            RG=0;
            RG=Radmodel.R(Srad(nr));
            RR= blkdiag(RR,RG);
            if paras.pseudoupdate==true
                ym=-1234.1234;
            else
                yjk=Radmodel.h(SatState{Ns}(Tk,:),Srad(nr))+sqrtm(Radmodel.R(Srad(nr)))*randn(Radmodel.hn,1);
                ym=vertcat(ym,yjk);
            end
            
        end
        
        
        
        
        if isempty(ym)==0 && sum(isnan(ym))==0
            [mz,Pz]=MeanCov(Y,w);
            Pz=Pz+RR;
            Pcc=CrossCov(x,mk,Y,mz,w);
            if paras.pseudoupdate==true
                K=Pcc/Pz;
                Pku=Pk-K*Pz*K';
                mku=mk;
                
            else
                disp(strcat('Measurement Update for sat ',num2str(Ns)))
                [mku,Pku]=KalmanUpdate(mk,Pk,mz,Pz,Pcc,ym);
            end
            
            
            SatState{Ns}.mu(Tk,:)=mku';
            SatState{Ns}.P(Tk,:)=reshape(Pku,1,36);
            
        else
            true
        end
        
    else      
        true
    end
    
    
end
