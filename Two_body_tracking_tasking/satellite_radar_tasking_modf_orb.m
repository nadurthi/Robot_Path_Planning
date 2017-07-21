clc
clear
close all
% everything ius done in the orbital spavce

Re=6378.1;
%% satellites
% Xsat0=getInitialrv_tles();
% Xsat0=getInitialrv_3D(20);
% Xsat0=getInitialrv_3D_3(20);
%  load('Initial20randOrbits_7')
%  load('Initialchinesedebris')
%  Xsat0([31,32,47],:)=[];
%  Xsat0=[Xsat0(1:15,:);Xsat1(1:30,:)];
%  Xsat0([ 10    11    12    13    14    15    17    18    21    22    23 25    28    29    33    34    37    41    43    44],:)=[];
%  
 

Xsat0_cart= [7200 0 0 1.0374090357 -1.3374090357 7.4771288355];
Xsat0= XYZ2OE(Xsat0_cart);

Nsat=size(Xsat0,1);
Radmodel.Nsat=Nsat;
%% radars (lat,long, altitude) in geod+edic
nsens_out=2;
close all
th=linspace(0,2*pi,5);
% phi=linspace(-pi/6,pi/6,5);
phi=0;
th(end)=[]; % perpendicular to equator
% phi(end)=[]; % is along the equator
% 0,0;;90,-60;270,60;270,-60
% Ang=[0,120;0,0;-60,-23;120,55;250,-10]*pi/180;
Ang=[-1.3800    1.9626
    0.1437   -0.0168
   -1.3865    1.8630
    0.1663   -3.0993
   -1.3936    1.5987
    0.1357   -0.0114
   -1.3910    1.7570
   -1.3879    1.3314
    0.1549   -3.0988
   -1.3842    1.9073
   -1.3818    1.2289];
 ind =[1 2 3 4  5 7 8 9 10 11]; % use 6 and 9
 Ang(ind,:)=[];
RadPos=zeros(size(Ang,1),3);
MeasCov=zeros(size(Ang,1),nsens_out^2);
k=1;
for i=1:1:size(Ang,1)
     RadPos(k,:)=[Ang(i,1),Ang(i,2),Re];  % [th,phi,Re] Re means on the surface
  R=blkdiag((0.5*pi/180)^2,(0.5*pi/180)^2);
%  R=blkdiag((0.1)^2);
  MeasCov(k,:)=reshape(R,1,nsens_out^2);
  k=k+1;
   
end
 Nrad=size(RadPos,1);
% p = randperm(25);
% p=p(1:10);
% MeasCov(p,1)=MeasCov(p,1)*100;
SensParas=[pi/6*ones(Nrad,1),2000*ones(Nrad,1)];  %misc paras such as [cone_angle,max dist of meas]
% SensParas(p,1)=pi/8;
% SensParas(p,2)=40000;
% % RadPos=[repmat(a,Nrad,2) + repmat((b-a),Nrad,2).*rand(Nrad,2),zeros(Nrad,1)];
% RadPos=[zeros(3,1),linspace(0,2*pi-pi/4,3)',Re*ones(3,1)];%[lat,long,r]
% RadPos=vertcat(RadPos,[linspace(0,2*pi-pi/6,10)',linspace(2*pi-pi/6,0,10)',Re*ones(10,1)]);
% % RadPos=vertcat(RadPos,[-pi/3*ones(5,1),linspace(pi/3,2*pi,5)',Re*ones(5,1)]);
% RadPos=vertcat(RadPos,[pi/2+pi/20,0,Re]);
% RadPos=vertcat(RadPos,[-pi/2-pi/20,0,Re]);

% RadRot_ecef2local=get_rad_rot_ecef2local(RadPos);

Radmodel.SensParas=SensParas;
Radmodel.RadPos=RadPos;
Radmodel.Nrad=Nrad;
Radmodel.hn=nsens_out;
% Radmodel.Q=blkdiag(0.01^2,0.01^2,0.01^2,1e-10,1e-10,1e-10);
Radmodel.Q=zeros(6,6);
Radmodel.G=@(x,Srad)radar_sens_penalty_orb(x,Srad,Radmodel,100);
Radmodel.h=@(x,Srad)radar_sens_orb(x,Srad,Radmodel);

Radmodel.R=@(Srad)reshape(MeasCov(Srad,:),Radmodel.hn,Radmodel.hn);

% plot_sat_radar_system2(yplottruth,Radmodel)
% xlabel('x')
% ylabel('y')
% zlabel('z')
% hold on
% plot3(yplottruth{14,1}(:,1),yplottruth{14,1}(:,2),yplottruth{14,1}(:,3),'ob--')
% hold off
%% time step setup
tf=1*24*60*60; % 48 hours 
dt=2*60; % 10 mins time step
t0=0;
Tvec=t0:dt:tf;
nt=length(t0:dt:tf);
MeasFreq=100; % i.e. after every 100 time steps (or dt)

%% generating the truth
% figure
% hold on

P0_cart=blkdiag(0.01,0.01,0.01,1e-8,1e-8,1e-8);
[X,w]=GH_points(Xsat0_cart',P0_cart,6);
for i=1:1:length(w)
   X(i,:)=XYZ2OE(X(i,:)); 
end
[mu0_orb,P0_orb]=MeanCov(X,w);

opt = odeset('reltol',1e-12,'abstol',1e-12);
ytruth=cell(Nsat,1);
yplottruth=cell(Nsat,1);

for i=1:Nsat

    [X,w]=UT_sigmapoints(Xsat0(i,:)',P0_orb,2);
    [a,b]=max(sqrt(sum(X.^2,2)));
    pprr=X(b,:);
[tt,xx]= twoBody_orbin(Tvec,pprr');
ytruth{i}=xx;
[tt,xx]=ode45(@twoBody,linspace(t0,tf,500),Xsat0_cart,opt);
yplottruth{i}=xx;
 i

end
ymeas=cell(Nsat,Nrad,nt);
for i=1:Nsat
    for j=1:1:Nrad
        for k=1:1:nt
ymeas{i,j,k}=Radmodel.h(ytruth{i,1}(k,:),j)+sqrtm(Radmodel.R(j))*randn(Radmodel.hn,1);
        end
    end
end


plot_sat_radar_system2(yplottruth,Radmodel)


% for k=1:1:nt
%     plot_sat_radar_system2(yplottruth,Radmodel)
%     hold on
%     for i=1:1:Nsat
%     X=ytruth{i}(k,:);
%     if norm(X(1:3))<Re
%         [k,i]
%     end
%     plot3(X(:,1),X(:,2),X(:,3),'bo')
%     
%     end
% %     axis([ -20000       15000      -12000        8000      -20000       20000])
%     pause(0.01)
% %     view(3,30)
%     hold off
%     clf
% end

%% checking if all the orbits are observable
Satobserve=zeros(Nsat,1);
for nsat=1:1:Nsat
    for nrad=1:1:Nrad
        for tt=1:1:length(Tvec)
            [yy,hh]=Radmodel.G(ytruth{nsat}(tt,:),nrad);
            if isnan(hh)==0
                Satobserve(nsat)=Satobserve(nsat)+1;
            end
        end
    end
end
Satobserve
%% Methods that are being used
% UT and CUT8
%generate initial sigma points for each sat

Xsat0;
XsigSat_ut=cell(Nsat,3);  %first is mean over time, second is cov over time
XsigSat_cut8=cell(Nsat,3);
H=zeros(1,Nsat);

SigSet_ut=cell(Nsat,3); %stores the latest points and weights so we donot have to ode45 from previous measurement update
SigSet_cut8=cell(Nsat,3);

for i=1:1:Nsat
    
    XsigSat_ut{i,1}=zeros(nt,6);
    XsigSat_cut8{i,1}=zeros(nt,6);
    XsigSat_ut{i,2}=zeros(nt,36);
    XsigSat_cut8{i,2}=zeros(nt,36);

    XsigSat_ut{i,3}=zeros(nt,1); % this is for meas time record
    XsigSat_cut8{i,3}=zeros(nt,1);% this is for meas time record
    
%     [X,w]=conjugate_dir_gausspts_till_8moment(Xsat0(i,:)',P0);
%     [a,b]=max(sqrt(sum(X.^2,2)));
%     pprr=X(b,:);
%     pprr=mvnrnd(Xsat0(i,:),P0);
    pprr=Xsat0(i,:);
    XsigSat_ut{i,1}(1,:)=pprr;
    XsigSat_cut8{i,1}(1,:)=pprr;
    XsigSat_ut{i,2}(1,:)=reshape(P0,1,36);
    XsigSat_cut8{i,2}(1,:)=reshape(P0,1,36);
    
    XsigSat_ut{i,3}(1)=1; 
    XsigSat_cut8{i,3}(1)=1;
    H(1,i)=log(det(P0));
    
    [X,w]=UT_sigmapoints(XsigSat_ut{i,1}(1,:)',P0,2);
    SigSet_ut{i,1}=X;
    SigSet_ut{i,2}=w;
    SigSet_ut{i,3}=1;
    tic
    [X,w]=conjugate_dir_gausspts_till_8moment(XsigSat_cut8{i,1}(1,:)',P0);
    SigSet_cut8{i,1}=X;
    SigSet_cut8{i,2}=w;
    SigSet_cut8{i,3}=1;% this point set corresponds to time step 1
    toc
end

%     
    
%% Propagate, task and filter
MeasPairs_ut=cell(nt,1);
MeasPairs_cut8=cell(nt,1);
MeasPairs_ut{1}=[[1:Nrad]',[1:Nrad]'];
MeasPairs_cut8{1}=[[1:Nrad]',[1:Nrad]'];






% figure(1)
% plot_sat_radar_system(XsigSat_ut,1,Radmodel,MeasPairs_ut,yplottruth)

% keyboard
for k=2:1:nt
    k
% figure(2) 
% clf
tic
 %Generate the MeasPairs
     [MeasPairs_ut,SigSet_ut]=SatRadarPair_modf_orb(MeasPairs_ut,XsigSat_ut,SigSet_ut,Radmodel,k,k,Tvec,'ut','MI');
     [MeasPairs_cut8,SigSet_cut8]=SatRadarPair_modf_orb(MeasPairs_cut8,XsigSat_cut8,SigSet_cut8,Radmodel,k,k,Tvec,'cut8','MI');
if rem(floor(k/MeasFreq),2)==0
for i=1:Nsat
    for j=1:1:Nrad
     ymeas{i,j,k}=ymeas{i,j,k}*NaN;
    end
end
end

% Measurement update only for the MeasPairs{k}
     [XsigSat_ut,SigSet_ut]=Meas_Update_mu_P_all_sats_modf2(XsigSat_ut,MeasPairs_ut,SigSet_ut,Radmodel,k,Tvec,ymeas,'ut',ytruth);
     [XsigSat_cut8,SigSet_cut8]=Meas_Update_mu_P_all_sats_modf2(XsigSat_cut8,MeasPairs_cut8,SigSet_cut8,Radmodel,k,Tvec,ymeas,'cut8',ytruth);
toc
% figure(2) 
% plot_sat_radar_system(XsigSat_ut,k,Radmodel,MeasPairs_ut,yplottruth)

% saveas(gca,strcat('SatRadpairUPDT_ut_',num2str(k)),'jpg')
% saveas(gca,strcat('SatRadpairUPDT_ut_',num2str(k)),'fig')

%      figure(3)    
% plot_sat_radar_system(XsigSat_cut8,k,Radmodel,MeasPairs_cut8,yplottruth)

% saveas(gca,strcat('SatRadpairUPDT_cut8_',num2str(k)),'jpg')
% saveas(gca,strcat('SatRadpairUPDT_cut8_',num2str(k)),'fig')
     
% figure(4)
% plot_enttropy_sats(XsigSat_ut,k,'r',Tvec);
% title('UT')
% figure(5)
% plot_enttropy_sats(XsigSat_cut8,k,'b',Tvec);
% title('CUT8')
% 
% figure(6)
Hut=plot_enttropy_sats_CUMULATE(XsigSat_ut,k,'r',Tvec);


% figure(7)
Hcut8=plot_enttropy_sats_CUMULATE(XsigSat_cut8,k,'b',Tvec);

figure(6)
plot(Tvec(1:k),Hcut8,Tvec(1:k),Hut,'--')
legend('cut','ut')

% figure(8)
Huterr=plot_est_err_sats_CUMULATE(XsigSat_ut,k,'r',Tvec,ytruth);


% figure(9)
Hcut8err=plot_est_err_sats_CUMULATE(XsigSat_cut8,k,'b',Tvec,ytruth);
figure(7)
plot(Tvec(1:k)/3600,Hcut8err(:,1),Tvec(1:k)/3600,Huterr(:,1),'--')
legend('cut','ut')

[Hcut8err(k,1),Huterr(k,1)]



    pause(0.05)
end

 save('Sat3Dtasking_estimation_300T','XsigSat_ut','XsigSat_cut8','MeasPairs_ut','MeasPairs_cut8','Hut','Hcut8','Huterr','Hcut8err','Tvec','ytruth','Xsat0','Radmodel','P0','MeasCov')