function Satellites=Propagate_sattask(Satellites,Constants,Tk,Tk1,method)

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


Nsats=Constants.Nsat;


opt = odeset('reltol',1e-12,'abstol',1e-12);

%%
parfor Ns=1:Nsats
    
    mk=Satellites{Ns}.mu(Tk,:)';
    Pk=reshape(Satellites{Ns}.P(Tk,:),Satellites{Ns}.fn,Satellites{Ns}.fn);
    
    
    [x,w]=qd_pts(mk,Pk);
    N=size(x,1);
    
    
    [N,n]=size(x);
    tT=Constants.Tvec(Tk:Tk1);
    Y=cell(N,1);
    F=Satellites{Ns}.f;
    for i=1:N
        [~,xx]=ode45(F ,tT,x(i,:)',opt);
        Y{i}=xx(end,:);
    end
    
    Z=zeros(size(x));
    for j=1:N
        Z(j,:)=Y{j};
    end
    
    W=repmat(w,1,n);
    mk1=sum(W.*Z,1)';
    
    
    MU=repmat(mk1',N,1);
    Z=Z-MU;
    Pk1=Z'*(W.*Z);
    
    Satellites{Ns}.mu(Tk1,:)=mk1;
    Satellites{Ns}.P(Tk1,:)=reshape(Pk1,1,Satellites{Ns}.fn*Satellites{Ns}.fn);
    
    
    
end


end