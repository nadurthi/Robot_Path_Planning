function [XsigSat,NewSigSet]=Meas_Update_mu_P_all_sats_modf_orb(XsigSat,MeasPairs,SigSet,Radmodel,k,Tvec,ymeas,method,ytruth)


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

    x=SigSet{Ns,1};
    w=SigSet{Ns,2};

[N,n]=size(x);
tT=Tvec([SigSet{Ns,3}:k]);
if k>SigSet{Ns,3}
x=propagate_sigma_pts_2body(x,tT);
end
[mk,Pk]=MeanCov(x,w);

xorb = XYZ2OE_multiple(x); % convert to orbital elements
[mkorb,Pkorb]=MeanCov(xorb,w);
%  Pk=Pk+Radmodel.Q;
%  [x,w]=qd_pts(mk,Pk);

    if sum(size(MeasPairs{k}))==0
        ind=[];
    else
    ind=find(MeasPairs{k}(:,1)==Ns);
    end
    
    if isempty(ind)==0
        Srad=MeasPairs{k}(ind,2);
%         Y=zeros(size(x,1),Radmodel.hn*length(Srad));
        Y=[];
        ym=[];
        RR=[];
        GG=[];
        for nr=1:1:length(Srad)
            flag1=0;
            ZZ=zeros(N,Radmodel.hn);
            G=zeros(N,1);
            H=zeros(N,1);
            for msi=1:1:N
                ZZ(msi,:)=Radmodel.h(x(msi,:)',Srad(nr));
                [gg,hh]=Radmodel.G(x(msi,:)',Srad(nr));
                G(msi)=gg;
                H(msi)=hh;
                %                 if isnan(ZZ(msi,1))==1
                %                     flag1=1;
                %                     break;
                %                 end
            end
            %             if flag1==0
            %             Y(msi,(nr-1)*Radmodel.hn+1:(nr*Radmodel.hn))=ZZ(:)';
            if sum(isnan(H))<length(H)/2
                Y=horzcat(Y,ZZ);
                RG=0;
                for ii=1:1:N
                    RG=RG+w(ii)*G(ii)^2*Radmodel.R(Srad(nr));
                end
                RR= blkdiag(RR,RG);
                ym=vertcat(ym,ymeas{Ns,Srad(nr),k});
                
            end
        end

        
        

        if isempty(ym)==0 && sum(isnan(ym))==0
        [mz,Pz]=MeanCov(Y,w);    
        Pz=Pz+RR;    
%         Pcc=CrossCov(x,mk,Y,mz,w);
        Pccorb=CrossCov(xorb,mkorb,Y,mz,w);
        [xuorb,Puorb]=KalmanUpdate(mkorb,Pkorb,mz,Pz,Pccorb,ym);
        [xu,Pu]=ORBmuCov_2_cartmuCov(xuorb,Puorb);
        XS1{Ns}(k,:)=xu(:)';
        XS2{Ns}(k,:)=reshape(Pu,1,length(xu)^2);
        XS3{Ns}(k)=1;
        
        [xk1orb,w]=qd_pts(xuorb,Puorb);
        x = OE2XYZ_multiple(xk1orb);
        
        NewSigSetpts{Ns,1}=x;
        NewSigSetwts{Ns,1}=w;
        disp('Meas Updated')
        else
            disp('NOOOO Meas Avoided')
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