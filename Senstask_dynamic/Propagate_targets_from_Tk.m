function Targets=Propagate_targets_from_Tk(Targets,time,Tk,method)
% propagate from Tk to Tk+1
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

%%
for Ntarg=1:Targets.Ntargs
    
    mk=Targets.xf{Ntarg}(Tk,:)';
    Pk=reshape(Targets.Pf{Ntarg}(Tk,:),Targets.fn(Ntarg),Targets.fn(Ntarg));
    
    if any(isreal(mk)==0) | any(isnan(mk)==1) | any(isreal(Pk)==0) | any(isnan(Pk)==1) | any(eig(Pk)<0)
        keyboard
    end
    
    [x,w]=qd_pts(mk,Pk);
    
    [N,n]=size(x);
    Z=zeros(size(x));
    F=Targets.f{Ntarg};
    for i=1:N
        Z(i,:)=F(time.dt,x(i,:)');
    end

    
    W=repmat(w,1,n);
    mk1=sum(W.*Z,1)';
    
    
    MU=repmat(mk1',N,1);
    Z=Z-MU;
    Pk1=Z'*(W.*Z)+Targets.Q{Ntarg};
    
    if any(isreal(mk1)==0) | any(isnan(mk1)==1) | any(isreal(Pk1)==0) | any(isnan(Pk1)==1) | any(eig(Pk1)<0)
        keyboard
    end
    
    
    Targets.xf{Ntarg}(Tk+1,:)=mk1;
    Targets.Pf{Ntarg}(Tk+1,:)=reshape(Pk1,1,Targets.fn(Ntarg)^2);
    
    
    
end


end