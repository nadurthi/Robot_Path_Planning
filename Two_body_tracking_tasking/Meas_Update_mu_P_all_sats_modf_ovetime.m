function [XsigSat,NewSigSet]=Meas_Update_mu_P_all_sats_modf_ovetime(XsigSat,MeasPairs,SigSet,Radmodel,k,Tvec,method,ytruth)
MeasPairs{k}
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
    
    XS1=cell(size(XsigSat,1),1);
    XS2=cell(size(XsigSat,1),1);
    XS3=cell(size(XsigSat,1),1);
    
for i=1:1:size(XsigSat,1)
    XS1{i}=XsigSat{i,1};
    XS2{i}=XsigSat{i,2};
    XS3{i}=XsigSat{i,3};
end
Nsatt=size(XsigSat,1);

NewSigSet=cell(size(SigSet,1),3);
NewSigSetpts=cell(size(SigSet,1),1);
NewSigSetwts=cell(size(SigSet,1),1);

for Ns=1:1:Nsatt
    

%         from_step=find(XS3{Ns}==1, 1, 'last' );

    
%     mu=XS1{Ns}(from_step,:)';
%     P=reshape(XS2{Ns}(from_step,:)',length(mu),length(mu));
    
%     [x,w]=qd_pts(mu,P);
    x=SigSet{Ns,1};
    w=SigSet{Ns,2};

[N,n]=size(x);
tT=Tvec([SigSet{Ns,3}:k]);
if k>SigSet{Ns,3}
% Y=cell(N,1);
    for sdd=1:1:N
    [tt,xx]=ode45(@twoBody,tT,x(sdd,:)',opt);
%     [tt,xx]=twoBodyKeplerProp(tT,x(sdd,:)');
    x(sdd,:)=xx(end,:);
    end
end

%     Z=zeros(size(x));
% %  for i=2:1:length(tT)
%         for j=1:1:N
%         Z(j,:)=x(j,:);
%         end
      
    W=repmat(w,1,n);
    mk=sum(W.*x,1)';
 

MU=repmat(mk',N,1);
xd=x-MU;
Pk=xd'*(W.*xd);

if max(size(MeasPairs{k}))==0
        XS1{Ns}(k,:)=mk(:)';
        XS2{Ns}(k,:)=reshape(Pk,1,length(mk)^2);
        NewSigSetpts{Ns,1}=x;
        NewSigSetwts{Ns,1}=w; 
    continue;
end
    ind=find(MeasPairs{k}(:,1)==Ns);
    
    if isempty(ind)==0
        Srad=MeasPairs{k}(ind,2);
%         [x,w]=qd_pts(mk,Pk);
%         [N,n]=size(x);
        Y=zeros(size(x,1),Radmodel.hn*length(Srad));
        flag1=0;
        for msi=1:1:N
            for nr=1:1:length(Srad)
                ZZ=Radmodel.h(x(msi,:)',Srad(nr));
                if isnan(ZZ(1))==1
                    flag1=1;
                end
                Y(msi,(nr-1)*Radmodel.hn+1:(nr*Radmodel.hn))=ZZ(:)';
            end
        end
        RR=[];
        for nr=1:1:length(Srad)
              RR= blkdiag(RR,Radmodel.R(Srad(nr)));
        end
            
        [N,nz]=size(Y);
        W=repmat(w,1,nz);
        mz=sum(W.*Y,1)';
 

        MU=repmat(mz',N,1);
        X=Y-MU;
        Pz=X'*(W.*X);
        
        Pz=Pz+RR;
        
        %cross cov
%         [x,w]=qd_pts(mk,Pk);
        Pcc=0;
        for msi=1:1:length(w)
            ZZ=[];
            for nr=1:1:length(Srad)
                zz=Radmodel.h(x(msi,:)',Srad(nr));
                ZZ=vertcat(ZZ,zz(:));
            end
            Pcc=Pcc+w(msi)*(x(msi,:)'-mk)*(ZZ-mz)';
        end
        %kalman gain
        K=Pcc/Pz;
        %update
        ym=[];
        Ytru=ytruth{Ns,1}(k,:);
        for nr=1:1:length(Srad)
            ym=vertcat(ym,Radmodel.h(Ytru',Srad(nr))+sqrtm(Radmodel.R(Srad(nr)))*randn(Radmodel.hn,1));
        end
        if isnan(sum(ym))>0
            flag1=1;
        end
%         figure(2)
%         plot3(Ytru(1),Ytru(2),Ytru(3),'ko',x(:,1),x(:,2),x(:,3),'ro')
%         legend('Ytrue','Post sigma pts')
%         grid
%         keyboard
        
        if flag1==0
        xu=mk+K*(ym-mz);
        Pu=Pk-K*Pz*K';
        XS1{Ns}(k,:)=xu(:)';
        XS2{Ns}(k,:)=reshape(Pu,1,length(xu)^2);
        XS3{Ns}(k)=1;
        [x,w]=qd_pts(xu,Pu);
        NewSigSetpts{Ns,1}=x;
        NewSigSetwts{Ns,1}=w;
        disp(strcat('updated Sat : ',num2str(Ns),' Cov = ',num2str(trace(Pu))));
        else
        XS1{Ns}(k,:)=mk(:)';
        XS2{Ns}(k,:)=reshape(Pk,1,length(mk)^2);
        NewSigSetpts{Ns,1}=x;
        NewSigSetwts{Ns,1}=w;    
        end
    else
        XS1{Ns}(k,:)=mk(:)';
        XS2{Ns}(k,:)=reshape(Pk,1,length(mk)^2);
        NewSigSetpts{Ns,1}=x;
        NewSigSetwts{Ns,1}=w;
    end
% end

end
for i=1:1:size(XsigSat,1)
  NewSigSet{i,1}=  NewSigSetpts{i};
  NewSigSet{i,2}=  NewSigSetwts{i};
  NewSigSet{i,3}=  k;
  
end

for i=1:1:size(XsigSat,1)
    XsigSat{i,1}=XS1{i};
    XsigSat{i,2}=XS2{i};
    XsigSat{i,3}=XS3{i};
end