function [MeasPairs,NewSigSet]=SatRadarPair_modf_ovetime(MeasPairs,XsigSat,SigSet,Radmodel,curk,fink,Tvec,method)



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



    XS1=cell(size(XsigSat,1),1);
    XS2=cell(size(XsigSat,1),1);
    XS3=cell(size(XsigSat,1),1);
    
for i=1:1:size(XsigSat,1)
    XS1{i}=XsigSat{i,1};
    XS2{i}=XsigSat{i,2};
    XS3{i}=XsigSat{i,3};
end

NewSigSet=cell(size(SigSet,1),3);
FIM0=cell(size(SigSet,1),1);
Muk=cell(size(SigSet,1),1);
invPhis=cell(size(SigSet,1),length(Tvec([curk:fink])));


for Ns=1:1:size(XsigSat,1)

    x=SigSet{Ns,1};
    w=SigSet{Ns,2};

    [N,n]=size(x);

%     tT=Tvec([SigSet{Ns,3},curk:fink]);
%     Y=cell(N,1);

if curk>SigSet{Ns,3}
    for i=1:1:N
%         mm=(tT(1)+tT(2))/2;
%         mm=[tT(1),mm,tT(2:end)];
%         [tt,xx]=twoBodyKeplerProp(mm,x(i,:)');
        [tt,xx]=ode45(@twoBody,Tvec([SigSet{Ns,3}:curk]),x(i,:)',opt);
        x(i,:)=xx(end,:);
%         if length(tT)==2
%         Y{i}=xx([1,end],:);
%         else

%         end
    end
end
NewSigSet{Ns,1}=x;

    W=repmat(w,1,n);
    mk=sum(W.*x,1)';
    MU=repmat(mk',N,1);
    x=x-MU;
    Pk=x'*(W.*x);
    FIM0{Ns}=inv(Pk);
    
    [tt,xx]=ode45(@twoBody,Tvec([curk:fink]),mk,opt);
    Muk{Ns}=xx;
    tT=Tvec([curk:fink]);
    phik_k1=STMcovPROP(mk,tT);
    for k=1:1:length(Tvec([curk:fink]))-1
        phi=reshape(phik_k1(k+1,:),6,6)*inv(reshape(phik_k1(k,:),6,6));
        invPhis{Ns,k}=inv(phi);
    end
    invPhis{Ns,length(Tvec([curk:fink]))}=eye(6);
end

tT=Tvec([curk:fink]);
nrads=Radmodel.Nrad;
nsats=size(SigSet,1);
nt=length(tT);


A=zeros(nt*nrads,nt*nrads*nsats);
H=cell(nsats,nrads,nt);
C=zeros(1,nt*nrads*nsats);
B=zeros(nsats,nrads,nt);
ss=1;
pp=1;
for i=1:1:nsats
    for k=1:1:nt
        for j=1:1:nrads
            H{i,j,k}=radar_sens_jac(Muk{i}(k,:),j,Radmodel);
            if isnan(radar_sens(Muk{i}(k,:),j,Radmodel))
                  C(pp,ss)=1; %this elemtn or mu variable has to be zero 
                 pp=pp+1;
            end
            B(i,j,k)=ss;
            ss=ss+1;
        end
    end
end

invR=cell(nrads,1);
 for j=1:1:nrads
            invR{j}=inv(Radmodel.R(j));
 end
 
 ss=1;       
 for k=1:1:nt
     for j=1:1:nrads
         for i=1:1:nsats
         A(ss,B(i,j,k))=1;
         end
         ss=ss+1;
     end
 end
b=ones(nt*nrads,1);
d=zeros(size(C,1),1);
% keyboard
cvx_begin 
    variable mmu(nt*nrads*nsats) binary
    maximize(FIMlinear_ovetime_cost(mmu,FIM0,invPhis,H,invR,nt,nrads,nsats))
    subject to
        A*mmu<=b
        C*mmu==d
cvx_end


TT=curk:fink;
for k=1:1:nt
    MeasPairs{TT(k)}=[];
end
ss=1;
for i=1:1:nsats
    for k=1:1:nt
        for j=1:1:nrads
            if mmu(ss)==1
                MeasPairs{TT(k)}=vertcat(MeasPairs{TT(k)},[i,j]);
            end
            ss=ss+1;
        end
    end
end

for i=1:1:size(XsigSat,1)
   NewSigSet{i,2}=SigSet{i,2};
   NewSigSet{i,3}=curk;
end


for k=1:1:nt
    SP=MeasPairs{TT(k)};
%     if length(unique(SP(:,2)))~=size(SP,1)
%         k
%     end
    for sp=1:1:size(SP,1)
        satt=SP(sp,1);
        rrad=SP(sp,2);
        if  isnan(radar_sens(Muk{satt}(k,:),rrad,Radmodel))
           MeasPairs{TT(k)}(sp,:)=[];
           disp('removing NaN assignment measurements')
        end
    end
    
end