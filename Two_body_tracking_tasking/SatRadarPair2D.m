function MeasPairs=SatRadarPair2D(MeasPairs,XsigSat,Radmodel,curk,fink,Tvec,method,FIMmethod)


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
FIMG=cell(size(XsigSat,1),1);


    XS1=cell(size(XsigSat,1),1);
    XS2=cell(size(XsigSat,1),1);
    XS3=cell(size(XsigSat,1),1);
    
for i=1:1:size(XsigSat,1)
    XS1{i}=XsigSat{i,1};
    XS2{i}=XsigSat{i,2};
    XS3{i}=XsigSat{i,3};
end

parfor Ns=1:1:size(XsigSat,1)
    FIMGinter=zeros(Radmodel.Nrad,length(curk:fink));
    
    from_step=find(XS3{Ns}==1, 1, 'last' );
    
    mu=XS1{Ns}(from_step,:)';
    P=reshape(XS2{Ns}(from_step,:)',length(mu),length(mu));
    
    [x,w]=qd_pts(mu,P);
    
    [N,n]=size(x);

    tT=Tvec([from_step,curk:fink]);
    Y=cell(N,1);
    for i=1:1:N
        mm=(tT(1)+tT(2))/2;
        mm=[tT(1),mm,tT(2:end)];
        [tt,xx]=ode45(@twoBody2D,mm,x(i,:)');
%         if length(tT)==2
%         Y{i}=xx([1,end],:);
%         else
        Y{i}=xx(3:end,:);
%         end
    end
%     keyboard
    tT=Tvec([curk:fink]);
    
     for tt=1:1:length(tT)
         Z=zeros(size(x));
         for j=1:1:N
             Z(j,:)=Y{j}(tt,:);
         end
         
         W=repmat(w,1,n);
         mk=sum(W.*Z,1)';
         
         MU=repmat(mk',N,1);
         Z=Z-MU;
         Pk=Z'*(W.*Z)+Radmodel.Q;
         
        [x,w]=qd_pts(mk,Pk);
        [N,n]=size(x);
        
         for Nr=1:1:Radmodel.Nrad
             
             Z=zeros(N,Radmodel.hn);
             for i=1:1:N
                 ZZ=Radmodel.h(x(i,:)',Nr);
                 if isnan(ZZ)
                     FIMGinter(Nr,tt)=-1e14;
                 else
                     Z(i,:)=ZZ(:)';
                 end
             end
             if FIMGinter(Nr,tt)==-1e14
                 continue;
             end
             
             
             [N,nz]=size(Z);
             W=repmat(w,1,nz);
             mz=sum(W.*Z,1)';
             
             MU=repmat(mz',N,1);
             Z=Z-MU;
             Pz=Z'*(W.*Z);
             
             Pz=Pz+Radmodel.R(Nr);
             
             %cross cov
             %         [x,w]=qd_pts(mk,Pk);
             Pcc=0;
             for i=1:1:N
                 Pcc=Pcc+w(i)*(x(i,:)'-mk)*(Radmodel.h(x(i,:)',Nr)-mz)';
             end
             %kalman gain
             K=Pcc/Pz;
             %update
             Pu=Pk-K*Pz*K';
             % FIM of reduction in covariance
             if strcmp(FIMmethod,'det')
%                  FIMG(Ns,Nr,tt)=det(inv(Pu)-inv(Pk));
             elseif  strcmp(FIMmethod,'trace')
                    
                    FIMGinter(Nr,tt)=-0.5*log(det(Pu)/det(Pk));
%                  FIMG(Ns,Nr,tt)=trace(inv(Pu)-inv(Pk));
             end
             
             
         end
     end
     FIMG{Ns}=FIMGinter;
end
for i=1:1:size(XsigSat,1)
    XsigSat{i,1}=XS1{i};
    XsigSat{i,2}=XS2{i};
    XsigSat{i,3}=XS3{i};
end
FIMGnew=zeros(size(XsigSat,1),Radmodel.Nrad,length(curk:fink));
tT=Tvec([curk:fink]);
for Ns=1:1:size(XsigSat,1)    
     for tt=1:1:length(tT)
         for Nr=1:1:Radmodel.Nrad
         FIMGnew(Ns,:,:)=FIMG{Ns};
         end
     end
end
FIMG=FIMGnew;
%% integer linear programming to find the combinations 
FIMG;
f=zeros(prod(size(FIMG)),1);
F=zeros(prod(size(FIMG)),2);
k=1;
II=eye(Radmodel.Nrad);
A=repmat(II,length(tT),size(XsigSat,1));
b=ones(size(A,1),1);
for tt=1:1:length(tT)
      for Ns=1:1:size(XsigSat,1);
          for Nr=1:1:Radmodel.Nrad
            f(k)=FIMG(Ns,Nr,tt);
            F(k,:)=[Ns,Nr];
            k=k+1;
          end
      end
end

x = bintprog(-f,A,b);
tind=curk:fink;
s1=1;
s2=Radmodel.Nrad*size(XsigSat,1);
for tt=1:1:length(tT)
    X=x(s1:s2);
    length(find(x==1))
MeasPairs{tind(tt)}=F(find(X==1),:);
s1=s2+1;
s2=s1+Radmodel.Nrad*size(XsigSat,1)-1;
end
F=FIMG(:,:,1);
end