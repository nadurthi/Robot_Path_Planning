function Satellites=MeasUpdate_sattask(MeasPairs,Satellites,Radars,Constants,Tk,ymeas,method,updatetype)

% MeasPairs{1}(i,j)
%               j=1,rad1  j=2,rad2   j=3,rad3  j=4,rad4   ...
% i=1 sat1
% i=2 sat2
% i=3 sat3
% .
% .
% .


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

RadarIds=1:1:Constants.Nrad;

parfor Ns=1:Constants.Nsat
    
    mk=Satellites{Ns}.mu(Tk,:)';
    Pk=reshape(Satellites{Ns}.P(Tk,:),Satellites{Ns}.fn,Satellites{Ns}.fn);
    
    
    
    RadarIds_forsat=MeasPairs{Tk}(Ns,:);
    RadarIds_forsat=RadarIds(RadarIds_forsat==1);
    
    
    
    if isempty(RadarIds_forsat)==0
        [x,w]=qd_pts(mk,Pk);
        N=size(x,1);
    
        %         Srad=MeasPairs{k}(ind,2);
        %         Y=zeros(size(x,1),Radmodel.hn*length(Srad));
        Y=[];
        ym=[];
        RR=[];
        GG=[];
        updateflag1=0;
        for nr=1:1:length(RadarIds_forsat)
            
            radid=RadarIds_forsat(nr);
            
            ZZ=zeros(N,Radars{radid}.hn);
            G=zeros(N,1);
            H=zeros(N,1);
            for msi=1:1:N
                ZZ(msi,:)=Radars{radid}.h( x(msi,:)' , Radars{radid}.PolarPositions, Radars{radid}.hn );
                [gg,hh]=Radars{radid}.G( x(msi,:)', Radars{radid}.PolarPositions, Radars{radid}.hn, Radars{radid}.ConeAngle,Radars{radid}.MaxRange,Radars{radid}.penalty);
                G(msi)=gg;
                H(msi)=hh;
                %                 if isnan(ZZ(msi,1))==1
                %                     flag1=1;
                %                     break;
                %                 end
            end
            %             if flag1==0
            %             Y(msi,(nr-1)*Radmodel.hn+1:(nr*Radmodel.hn))=ZZ(:)';
            
            if strcmp(updatetype,'pseudoupdate')==1
                % do pseudo update here
                
                Y=horzcat(Y,ZZ);
                RG=Radars{radid}.R;
                RR= blkdiag(RR,RG);
                
                updateflag1=1;

                
            else
                if sum(isnan(H))<length(H)/2
                    Y=horzcat(Y,ZZ);
                    %  Do not apply penalty
                    %                     RG=0;
                    %                     for ii=1:1:N
                    %                         RG=RG+w(ii)*G(ii)^2*Radars{radid}.R;
                    %                     end
                    RG=Radars{radid}.R;
                    RR= blkdiag(RR,RG);
<<<<<<< HEAD
                    ym=vertcat(ym,ymeas{Ns}{radid,Tk});
=======
                    ym=vertcat(ym,ymeas{Ns,radid,Tk});
>>>>>>> 30c240a9f0bda66371f3d696fac22e6427c3e422
                    updateflag1=1;

                end
                
            end
        end
        
        
        
        
        if updateflag1==1
            [mz,Pz]=MeanCov(Y,w);
            Pz=Pz+RR;
            Pcc=CrossCov(x,mk,Y,mz,w);
            if strcmp(updatetype,'pseudoupdate')==1
                [xu,Pu]=KalmanUpdate(mk,Pk,mz,Pz,Pcc,mz);
            else
                [xu,Pu]=KalmanUpdate(mk,Pk,mz,Pz,Pcc,ym);
            end
            
            Satellites{Ns}.mu(Tk,:)=xu;
            Satellites{Ns}.P(Tk,:)=reshape(Pu,1,Satellites{Ns}.fn*Satellites{Ns}.fn);
            
%             disp(strcat('MeasUpdated for satellite         :',' ',num2str(Ns) ) )
        else
%             disp(strcat('NOOOO Meas for satellite          :',' ',num2str(Ns) ) )
            
            
        end
    else
%         disp(strcat('NOOOO Radars Tasked for satellite :',' ',num2str(Ns) )  )
    end
    % end
    
end
