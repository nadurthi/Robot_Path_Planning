function R=check_rad_sat_pair_visibility(Nsat,Nrad,XsigSat,Radmodel,Tk,Tvec,method)


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


Tsim=Tvec(Tk);



for Ns=Nsat
    
    mk=XsigSat{Ns,1}(Tk,:)';
    Pk=reshape(XsigSat{Ns,2}(Tk,:),6,6);
    
    
    [x,w]=qd_pts(mk,Pk);
    N=size(x,1);
    
    
    
    Srad=Nrad;
    
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
        
        if sum(isnan(H))<length(H)/2
            R=1;
        else
            R=-1;
        end
    end

end
