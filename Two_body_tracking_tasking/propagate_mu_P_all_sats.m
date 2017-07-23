function XsigSat=propagate_mu_P_all_sats(XsigSat,from_step,to_step,Tvec,method)

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
for Ns=1:1:size(XsigSat,1);

    if strcmp(from_step,'last_meas_step')
        from_step=find(XsigSat{Ns,3}==1, 1, 'last' );
    end
    
    mu=XsigSat{Ns,1}(from_step,:)';
    P=reshape(XsigSat{Ns,2}(from_step,:)',length(mu),length(mu));
    
    [x,w]=qd_pts(mu,P);

[N,n]=size(x);
tT=Tvec([from_step:to_step]);
Y=cell(N,1);
    for i=1:1:N
    [tt,xx]=ode45(@twoBody,tT,x(i,:)',opt);
    if length(tT)==2
        Y{i}=xx([1,end],:);
    else
    Y{i}=xx;
    end
    end
    Z=zeros(size(x));
    for i=2:1:length(tT)
        for j=1:1:N
        Z(j,:)=Y{j}(i,:);
        end
      
    W=repmat(w,1,n);
    mk=sum(W.*Z,1)';
 

MU=repmat(mk',N,1);
Z=Z-MU;
Pk=Z'*(W.*Z);
XsigSat{Ns,1}(from_step+i-1,:)=mk(:)';
XsigSat{Ns,2}(from_step+i-1,:)=reshape(Pk,1,length(mu)^2);
    end
%     Ns
end
