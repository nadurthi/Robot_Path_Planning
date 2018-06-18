function Targets=MeasUpdate_targets(Targets,Sensors,Tk,method)
% for each target and do the measureemnts
% if there are no measurements then no measurement update



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


for Ntarg=1:Targets.Ntargs
    
    mk=Targets.xf{Ntarg}(Tk,:)';
    Pk=reshape(Targets.Pf{Ntarg}(Tk,:),Targets.fn(Ntarg),Targets.fn(Ntarg));
    
    if any(isreal(mk)==0) | any(isnan(mk)==1) | any(isreal(Pk)==0) | any(isnan(Pk)==1) | any(eig(Pk)<0)
        keyboard
    end
    
    
    
    [x,w]=qd_pts(mk,Pk);
    N=size(x,1);
    
    Y=[];
    ym=[];
    RR=[];
    updateflag1=0;
    sensors_used = zeros(1,Sensors.Nsens);
    YM=cell(1,Sensors.Nsens);
    for nr= 1:Sensors.Nsens
        [y,g,~] = Sensors.h{nr}( Targets.truth{Ntarg}(Tk,:) , Sensors.xf{nr}(Tk,:), Sensors.FOV{nr} );
        YM{nr}=y;
        if g==1
            sensors_used(nr)=1;
        end
    end
    
    est_sensors_used = zeros(1,Sensors.Nsens);
    ZZ=cell(1,Sensors.Nsens);
    for nr= 1:Sensors.Nsens
        
        ZZ{nr}=zeros(N,Sensors.hn(nr));
        G=zeros(N,1);
        for msi=1:1:N
            [y,g,~] = Sensors.h{nr}( x(msi,:)' , Sensors.xf{nr}(Tk,:), Sensors.FOV{nr} );
            ZZ{nr}(msi,:)=y;
            G(msi)=g;
        end
        
        
        if sum(G==-1)<length(G)/2
            est_sensors_used(nr)=1;
        end
    end

    sens_used_flags = sensors_used & est_sensors_used;
    
    if sum(sens_used_flags)==0
        updateflag1 = 0;
    else
        for nr= 1:Sensors.Nsens
            if sens_used_flags(nr)==1
                disp(['meas updated ',num2str(Ntarg),' ',num2str(nr),' ',num2str(Tk)])
                ym=vertcat(ym,YM{nr});
                Y=horzcat(Y,ZZ{nr});
                RR=blkdiag(RR,Sensors.R{nr});
            end
        end
        updateflag1 = 1;
    end
    
    
    
    if updateflag1==1
        [mz,Pz]=MeanCov(Y,w);
        Pz=Pz+RR;
        Pcc=CrossCov(x,mk,Y,mz,w);

        [xu,Pu]=KalmanUpdate(mk,Pk,mz,Pz,Pcc,ym);

        if any(isreal(xu)==0) | any(isnan(xu)==1) | any(isreal(Pu)==0) | any(isnan(Pu)==1) | any(eig(Pu)<0)
            keyboard
        end
        
        Targets.xf{Ntarg}(Tk,:)=xu;
        Targets.Pf{Ntarg}(Tk,:)=reshape(Pu,1,Targets.fn(Ntarg)^2);

        
    end
    

end


end
