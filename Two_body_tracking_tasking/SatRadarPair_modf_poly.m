function [MeasPairs,NewSigSet]=SatRadarPair_modf_poly(MeasPairs,XsigSat,SigSet,Radmodel,curk,fink,Tvec,method,FIMmethod)



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

Vis=zeros(Radmodel.Nsat,Radmodel.Nrad);
PX=cell(Radmodel.Nsat,1);
PXZ=cell(Radmodel.Nsat,1);
PZ=cell(Radmodel.Nsat,1);

for Ns=1:1:size(XsigSat,1)
    
    x=SigSet{Ns,1};
    w=SigSet{Ns,2};
    
    
    [N,n]=size(x);
    
    tT=Tvec([SigSet{Ns,3},curk]);
    
    if curk>SigSet{Ns,3}
        x=propagate_sigma_pts_2body(x,tT);
    end
    
    
    [mk,Pk]=MeanCov(x,w);
    %  Pk=Pk+Radmodel.Q;
    NewSigSet{Ns,1}=x;
    PX{Ns}=Pk;
    %   [x,w]=qd_pts(mk,Pk);
    %  NewSigSet{Ns,1}=x;
    %   NewSigSet{Ns,2}=w;
    Zpts=[];
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
        Zpts=horzcat(Zpts,Z);
        if sum(isnan(H))>length(H)/2
            Vis(Ns,Nr)=0;
        else
            Vis(Ns,Nr)=1;
        end
    end
    
    if sum(Vis(Ns,:))>0
        
        RR=[];
        for Nr=1:1:Radmodel.Nrad
            if Vis(Ns,Nr)==1
                RR=blkdiag(RR,Radmodel.R(Nr));
            else
                RR=blkdiag(RR,Radmodel.R(Nr));
            end
        end
        %         keyboard
        [mz,Pz]=MeanCov(Zpts,w);
        PZ{Ns}=Pz+RR;
        PXZ{Ns}=CrossCov(x,mk,Zpts,mz,w);
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

if sum(sum(Vis))==0
    return
end

Px=cell(1);
Pz=cell(1);
Pxz=cell(1);
k=1;
INDmat=[];
for Ns=1:1:size(XsigSat,1)
    
    if sum(Vis(Ns,:))>0
        Px{k}=PX{Ns}*1;
        Pz{k}=PZ{Ns}*1;
        Pxz{k}=PXZ{Ns}*1;
        INDmat=vertcat(INDmat,[k,Ns]);
        k=k+1;
    end
    
end
% keyboard
try
    mset clear;
    nx=length(Px);
    nz=Radmodel.Nrad;
    mpol('v',nx,nz);
    [p,q]=MI_num_dem_poly_zver_multi_x(v,Px,Pxz,Pz,Radmodel.hn,nz,nx,0);
    
    K=(ones(nx,nz)-v).*v;
    K=reshape(K,nx*nz,1);
    mpol('S',nz,1);
    mpol RS;
    RS=0;
    for i=1:1:nz
        S(i)=sum(v(:,i));
        RS=RS+S(i);
    end
    K=[K==0;S<=1;RS>=0.1];
    VM=Vis(INDmat(:,2),:);
    for i=1:1:size(VM,1)
        for j=1:1:size(VM,2)
            if  VM(i,j)==0
                K=[K;v(i,j)==0];
            end
        end
    end
    Prob = msdp(max(q),K);
    [status,obj] = msol(Prob) ;
    vsolz = round(double(v));
catch
    keyboard
    q=q/min(coef(q));
    list = listvar(q);
    y=coef(q);
    nc=length(y);
    [ys,ind]=sort(y);
    A=diag(-1*ones(1,nc))+diag(1*ones(1,nc-1),1);
    
    cvx_begin
    variable xx(nc) nonnegative
    minimize(norm(ys-xx,4)+0.9*norm(A*xx,2))
    subject to
    ys>=xx
    cvx_end
    
    yr=zeros(nc,1);
    for i=1:1:nc
        yr(i)=xx(ind==i) ;
    end
[y,yr]

    q=sum(list.^2.*yr);
    Prob = msdp(max(q),K);
    [status,obj] = msol(Prob) ;
    vsolz = round(double(v));
    
    
    
end
pp=1;
for k=1:1:nx
    for Nr=1:1:nz
        if vsolz(k,Nr)==1
            MeasPairs{curk}(pp,:)=[INDmat(k,2),Nr];
            pp=pp+1;
        end
    end
end

end