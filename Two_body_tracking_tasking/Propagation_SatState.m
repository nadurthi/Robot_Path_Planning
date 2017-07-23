function SatState=Propagation_SatState(Nsats,SatState,T0,Tk,Tvec,method,paras)
% Use T0 initial uncertainty to propogate upto toime Tk
% also load the intermediate means and covs to SatState

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

% Nsats=lenght(SatState);

Tsim=Tvec(T0:Tk);

% Mu_t=cell(Nsatt,1);
% P_t=cell(Nsatt,1);

opt = odeset('reltol',1e-12,'abstol',1e-12);
for Ns=Nsats
    mu= SatState{Ns}.mu(T0,:)';
    P=reshape(SatState{Ns}.P(T0,:),6,6);
    
    [X0,w]=qd_pts(mu,P);
    Xt=zeros(length(Tsim),6,length(w));
% keyboard
    parfor ni=1:length(w)
        
        if length(Tsim)==2
            [~,xx]=ode45(@twoBody,Tsim,X0(ni,:)',opt);
            xx=xx([1,end],:);
        else
        [~,xx]=ode45(@twoBody,Tsim,X0(ni,:)',opt);
        end
        Xt(:,:,ni)=xx;
    end

    % Now calculating allthe moments for each of the time steps
%     Mu_t{Ns}=zeros(length(Tsim),6);
%     P_t{Ns}=zeros(length(Tsim),36);
    for tt=2:1:length(Tsim)
        XX=zeros(size(Xt,3),size(Xt,2));
        for a=1:1:size(Xt,3)
            XX(a,:)=Xt(tt,:,a);
        end
%         keyboard
        [mk,Pk]=MeanCov(XX,w);
        if paras.ProcessNoise==true
            Pk=Pk+SatState.Q;
        else
            Pk=Pk;
        end
        
        SatState{Ns}.mu(T0+tt-1,:)=mk';
        SatState{Ns}.P(T0+tt-1,:)=reshape(Pk,1,36);
        
    end
    
end

