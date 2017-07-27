function [MeasPairs,NewSigSet]=SatRadarPair_modf(MeasPairs,XsigSat,SigSet,Radmodel,curk,fink,Tvec,method,FIMmethod)



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
% FIMG=zeros(size(XsigSat,1),Radmodel.Nrad,length(curk:fink));
FIMG=zeros(size(XsigSat,1),Radmodel.Nrad);


    XS1=cell(size(XsigSat,1),1);
    XS2=cell(size(XsigSat,1),1);
    XS3=cell(size(XsigSat,1),1);
    
for i=1:1:size(XsigSat,1)
    XS1{i}=XsigSat{i,1};
    XS2{i}=XsigSat{i,2};
    XS3{i}=XsigSat{i,3};
end
NewSigSet=cell(size(SigSet,1),3);

for Ns=1:1:size(XsigSat,1)

    x=SigSet{Ns,1};
    w=SigSet{Ns,2};

    
    [N,n]=size(x);

    tT=Tvec([SigSet{Ns,3},curk]);
    
if curk>SigSet{Ns,3}
x=propagate_sigma_pts_2body(x,tT);
end
NewSigSet{Ns,1}=x;

[mk,Pk]=MeanCov(x,w);
         
         for Nr=1:1:Radmodel.Nrad
             
             Z=zeros(N,Radmodel.hn);
             G=zeros(N,1);
             H=zeros(N,1);
             for i=1:1:N
                 ZZ=Radmodel.h(x(i,:)',Nr);
                 [gg,hh]=Radmodel.G(x(i,:)',Nr);
                 G(i)=gg;
                 H(i)=hh;
%                  if isnan(hh)==1
%                      FIMG(Ns,Nr)=-1e14;
%                      break;
%                  else
                     Z(i,:)=ZZ(:)';
%                  end
             end
             
             if sum(isnan(H))>length(H)/2
                 FIMG(Ns,Nr)=-1e14;
                 continue;
             end
             
             [mz,Pz]=MeanCov(Z,w);
             
             RR=0;
             for ii=1:1:N
                 RR=RR+w(ii)*G(ii)^2*Radmodel.R(Nr);
             end
             Pz=Pz+RR;

             Pcc=CrossCov(x,mk,Z,mz,w);

             [~,Pu]=KalmanUpdate(mk,Pk,mz,Pz,Pcc,NaN);
             
             % FIM of reduction in covariance
             if strcmp(FIMmethod,'det')
%                  FIMG(Ns,Nr,tt)=det(inv(Pu)-inv(Pk));
             elseif  strcmp(FIMmethod,'MI')
                     FIMG(Ns,Nr)=-0.5*log(det(Pu)/det(Pk));
             elseif  strcmp(FIMmethod,'prior')
                     FIMG(Ns,Nr)=trace(Pk);        
             elseif  strcmp(FIMmethod,'trace')
                    FIMG(Ns,Nr)=trace(inv(Pu));
%                     FIMGinter(Nr,tt)=trace(Pk);
%                     FIMGinter(Nr,tt)=-0.5*log(det(Pu)/det(Pk));
%                  FIMG(Ns,Nr,tt)=trace(inv(Pu)-inv(Pk));
             end
             
             
         end
     end


for i=1:1:size(XsigSat,1)
   NewSigSet{i,2}=SigSet{i,2};
   NewSigSet{i,3}=curk;
end


for i=1:1:size(XsigSat,1)
    XsigSat{i,1}=XS1{i};
    XsigSat{i,2}=XS2{i};
    XsigSat{i,3}=XS3{i};
end

%% integer linear programming to find the combinations 

f=zeros(prod(size(FIMG)),1);
F=zeros(prod(size(FIMG)),2);

II=eye(Radmodel.Nrad);
A=repmat(II,1,size(XsigSat,1));
b=ones(size(A,1),1);
k=1;
      for Ns=1:1:size(XsigSat,1);
          for Nr=1:1:Radmodel.Nrad
            f(k)=FIMG(Ns,Nr);
            F(k,:)=[Ns,Nr];
            k=k+1;
          end
      end

if min(-f)==1e14
    return
end

mmu = bintprog(-f,A,b);

MeasPairs{curk}=[];
for tt=1:1:prod(size(FIMG))
    if mmu(tt)==1
    MeasPairs{curk}=vertcat(MeasPairs{curk},F(tt,:));
    end
end

end