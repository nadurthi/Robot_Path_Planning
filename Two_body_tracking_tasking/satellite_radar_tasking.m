clc
clear
close all

Re=6378.1;
%% satellites
% Xsat0=getInitialrv_tles();
% Xsat0=getInitialrv_3D(20);
load('Initial20randOrbits')
% Xsat0([2,6],:)=[];
% Xsat0=Xsat0([1:1:10],:);
Nsat=size(Xsat0,1);
Radmodel.Nsat=Nsat;
%% radars (lat,long, altitude) in geod+edic
nsens_out=3;

th=linspace(-pi/2,pi/2,6);
phi=linspace(0,2*pi,6);
th(end)=[]; % perpendicular to equator
phi(end)=[]; % is along the equator

RadPos=zeros(25,3);
MeasCov=zeros(25,nsens_out^2);
k=1;
for i=1:1:5
    for j=1:1:5
  RadPos(k,:)=[th(i),phi(j),Re];  % [th,phi,Re] Re means on the surface
  R=blkdiag((0.1)^2,(1*pi/180)^2,(1*pi/180)^2);
  MeasCov(k,:)=reshape(R,1,nsens_out^2);
  k=k+1;
    end
end
Nrad=size(RadPos,1);
p = randperm(25);
p=p(1:10);
MeasCov(p,1)=MeasCov(p,1)*10^8;
SensParas=[pi/3*ones(Nrad,1),10000*ones(Nrad,1)];  %misc paras such as [cone_angle,max dist of meas]
SensParas(p,1)=pi/8;
SensParas(p,2)=40000;
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
Radmodel.Q=zeros(6);
Radmodel.h=@(x,Srad)radar_sens(x,Srad,Radmodel);

Radmodel.R=@(Srad)reshape(MeasCov(Srad,:),nsens_out,nsens_out);



%% time step setup
tf=48*60*60; % 48 hours 
dt=1*60; % 10 mins time step
t0=0;
Tvec=t0:dt:tf;
nt=length(t0:dt:tf);

%% generating the truth
% figure
% hold on
opt = odeset('reltol',1e-12,'abstol',1e-12);
ytruth=cell(Nsat,1);
yplottruth=cell(Nsat,1);
for i=1:Nsat
[tt,xx]=ode45(@twoBody,Tvec,Xsat0(i,:)',opt);
ytruth{i}=xx;
[tt,xx]=ode45(@twoBody,linspace(t0,tf,500),Xsat0(i,:)',opt);
yplottruth{i}=xx;
% i
% plot3(xx(:,1),xx(:,2),xx(:,3))
% keyboard
end

%% checking if all the orbits are observable
% Satobserve=zeros(Nsat,1);
% for nsat=1:1:Nsat
%     for nrad=1:1:Nrad
%         for tt=1:1:length(Tvec)
%             yy=Radmodel.h(ytruth{nsat}(tt,:),nrad);
%             if isnan(yy(1))==0
%                 Satobserve(nsat)=Satobserve(nsat)+1;
%             end
%         end
%     end
% end
% Satobserve
%% Methods that are being used
% UT and CUT8
%generate initial sigma points for each sat
P0=blkdiag(0.01,0.01,0.01,1e-09,1e-09,1e-09);
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
    
    XsigSat_ut{i,1}(1,:)=mvnrnd(Xsat0(i,:),2*P0);
    XsigSat_cut8{i,1}(1,:)=mvnrnd(Xsat0(i,:),2*P0);
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
     [MeasPairs_ut,SigSet_ut]=SatRadarPair_modf(MeasPairs_ut,XsigSat_ut,SigSet_ut,Radmodel,k,k,Tvec,'ut','MI');
     [MeasPairs_cut8,SigSet_cut8]=SatRadarPair_modf(MeasPairs_cut8,XsigSat_cut8,SigSet_cut8,Radmodel,k,k,Tvec,'cut8','MI');

% Measurement update only for the MeasPairs{k}
     [XsigSat_ut,SigSet_ut]=Meas_Update_mu_P_all_sats_modf(XsigSat_ut,MeasPairs_ut,SigSet_ut,Radmodel,k,Tvec,'ut',ytruth);
     [XsigSat_cut8,SigSet_cut8]=Meas_Update_mu_P_all_sats_modf(XsigSat_cut8,MeasPairs_cut8,SigSet_cut8,Radmodel,k,Tvec,'cut8',ytruth);
toc
% figure(2) 
% plot_sat_radar_system(XsigSat_ut,k,Radmodel,MeasPairs_ut,yplottruth)
% 
% saveas(gca,strcat('SatRadpairUPDT_ut_',num2str(k)),'jpg')
% saveas(gca,strcat('SatRadpairUPDT_ut_',num2str(k)),'fig')

%      figure(3)    
% plot_sat_radar_system(XsigSat_cut8,k,Radmodel,MeasPairs_cut8,yplottruth)

% saveas(gca,strcat('SatRadpairUPDT_cut8_',num2str(k)),'jpg')
% saveas(gca,strcat('SatRadpairUPDT_cut8_',num2str(k)),'fig')
     
figure(4)
plot_enttropy_sats(XsigSat_ut,k,'r',Tvec);
title('UT')
figure(5)
plot_enttropy_sats(XsigSat_cut8,k,'b',Tvec);
title('CUT8')
% 
figure(6)
Hut=plot_enttropy_sats_CUMULATE(XsigSat_ut,k,'r',Tvec);


figure(7)
Hcut8=plot_enttropy_sats_CUMULATE(XsigSat_cut8,k,'b',Tvec);

figure(8)
Huterr=plot_est_err_sats_CUMULATE(XsigSat_ut,k,'r',Tvec,ytruth);


figure(9)
Hcut8err=plot_est_err_sats_CUMULATE(XsigSat_cut8,k,'b',Tvec,ytruth);

     
    pause(0.05)
end

save('Sat3Dtasking_estimation_sim1','XsigSat_ut','XsigSat_cut8','MeasPairs_ut','MeasPairs_cut8','Hut','Hcut8')