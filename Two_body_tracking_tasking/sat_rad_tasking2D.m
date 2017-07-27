clc
clear
close all

Re=6378.1;
%% satellites
Xsat0=getInitialrv_2D(10);
% Xsat0([2,6],:)=[];
%  Xsat0=Xsat0([1:1:5],:);
Nsat=size(Xsat0,1);
Radmodel.Nsat=Nsat;
%% radars (lat,long, altitude) in geod+edic

nsens_out=2;  % i.e. range, bearing(from north), azimuth(from up)

RadPos=[];
% RadPos=[repmat(a,Nrad,2) + repmat((b-a),Nrad,2).*rand(Nrad,2),zeros(Nrad,1)];
RadPos=[linspace(0,2*pi-pi/4,5)',Re*ones(5,1)];%[lat,long,r]
% RadPos=vertcat(RadPos,[linspace(0,2*pi-pi/6,10)',linspace(2*pi-pi/6,0,10)',Re*ones(10,1)]);
% RadPos=vertcat(RadPos,[-pi/3*ones(5,1),linspace(pi/3,2*pi,5)',Re*ones(5,1)]);
% RadPos=vertcat(RadPos,[pi/2+pi/20,0,Re]);
% RadPos=vertcat(RadPos,[-pi/2-pi/20,0,Re]);
Nrad=size(RadPos,1);

% RadRot_ecef2local=get_rad_rot_ecef2local(RadPos);
SensParas=[[10,20,30,40,50,60]'*pi/180,linspace(40000,2000,6)'];  %misc paras such as [cone_angle,max dist of meas]
Radmodel.SensParas=SensParas;
Radmodel.RadPos=RadPos;
Radmodel.Nrad=Nrad;
Radmodel.hn=nsens_out;

Radmodel.h=@(x,Srad)radar_sens2D(x,Srad,Radmodel);

%covariances for meas
MeasCov=zeros(Nrad,nsens_out^2);
R=blkdiag((0.1)^2,(1*pi/180)^2);
for i=1:1:Nrad
    MeasCov(i,:)=reshape(i*R,1,nsens_out^2);
end

Radmodel.R=@(Srad)reshape(MeasCov(Srad,:),nsens_out,nsens_out);
Radmodel.Q=0;%blkdiag((0.00000005)^2,(0.00000005)^2,1e-20,1e-20);
%% time step setup
tf=48*60*60; % 48 hours 
dt=5*60; % 10 mins time step
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
[tt,xx]=ode45(@twoBody2D,Tvec,Xsat0(i,:)',opt);
ytruth{i}=xx;
[tt,xx]=ode45(@twoBody2D,linspace(t0,tf,500),Xsat0(i,:)',opt);
yplottruth{i}=xx;
% i
% plot3(xx(:,1),xx(:,2),xx(:,3))
% keyboard
end
%% Methods that are being used
% UT and CUT8
%generate initial sigma points for each sat
P0=blkdiag(0.01,0.01,1e-07,1e-07);
Xsat0;
XsigSat_ut=cell(Nsat,3);  %first is mean over time, second is cov over time
XsigSat_cut8=cell(Nsat,3);


for i=1:1:Nsat
    
    XsigSat_ut{i,1}=zeros(nt,4);
    XsigSat_cut8{i,1}=zeros(nt,4);
    XsigSat_ut{i,2}=zeros(nt,16);
    XsigSat_cut8{i,2}=zeros(nt,16);

    XsigSat_ut{i,3}=zeros(nt,1); % this is for meas time record
    XsigSat_cut8{i,3}=zeros(nt,1);% this is for meas time record
    
    XsigSat_ut{i,1}(1,:)=mvnrnd(Xsat0(i,:),P0);
    XsigSat_cut8{i,1}(1,:)=mvnrnd(Xsat0(i,:),P0);
    XsigSat_ut{i,2}(1,:)=reshape(P0,1,16);
    XsigSat_cut8{i,2}(1,:)=reshape(P0,1,16);
    
    XsigSat_ut{i,3}(1)=1; 
    XsigSat_cut8{i,3}(1)=1;

end

%     
    
%% Propagate, task and filter
MeasPairs_ut=cell(nt,1);
MeasPairs_cut8=cell(nt,1);
MeasPairs_ut{1}=[[1:Nrad]',[1:Nrad]'];
MeasPairs_cut8{1}=[[1:Nrad]',[1:Nrad]'];



figure(1)
plot_sat_radar_system2D(XsigSat_ut,1,Radmodel,MeasPairs_ut,yplottruth)

keyboard

for k=2:1:nt
    k
%   keyboard
    figure(2)
    clf
%     figure(3)
%     clf
    % Propagate for 1 step
%      XsigSat_ut=propagate_mu_P_all_sats2D(XsigSat_ut,'last_meas_step',k,Tvec,'ut',Radmodel);
%      XsigSat_cut8=propagate_mu_P_all_sats2D(XsigSat_cut8,'last_meas_step',k,Tvec,'cut8',Radmodel);
    
     %Generate the MeasPairs
     tic
     MeasPairs_ut=SatRadarPair2D(MeasPairs_ut,XsigSat_ut,Radmodel,k,k,Tvec,'ut','trace');
     MeasPairs_cut8=SatRadarPair2D(MeasPairs_cut8,XsigSat_cut8,Radmodel,k,k,Tvec,'cut8','trace');
     toc
     % Measurement update only for the MeasPairs{k}
     XsigSat_ut=Meas_Update_mu_P_all_sats2D(XsigSat_ut,MeasPairs_ut,Radmodel,k,Tvec,'ut',ytruth);
     XsigSat_cut8=Meas_Update_mu_P_all_sats2D(XsigSat_cut8,MeasPairs_cut8,Radmodel,k,Tvec,'cut8',ytruth);

%      figure(2)    
% plot_sat_radar_system2D(XsigSat_ut,k,Radmodel,MeasPairs_ut,yplottruth)
% 
% saveas(gca,strcat('SatRadpairUPDT_ut_',num2str(k)),'jpg')
% saveas(gca,strcat('SatRadpairUPDT_ut_',num2str(k)),'fig')
% 
%      figure(3)    
% plot_sat_radar_system2D(XsigSat_cut8,k,Radmodel,MeasPairs_cut8,yplottruth)
% 
% saveas(gca,strcat('SatRadpairUPDT_cut8_',num2str(k)),'jpg')
% saveas(gca,strcat('SatRadpairUPDT_cut8_',num2str(k)),'fig')
%      
% figure(4)
% plot_enttropy_sats(XsigSat_ut,k,'r',Tvec);
% title('UT')
% figure(5)
% plot_enttropy_sats(XsigSat_cut8,k,'b',Tvec);
% title('CUT8')
% 
% figure(6)
% Hut=plot_enttropy_sats_CUMULATE(XsigSat_ut,k,'r',Tvec);
% 
% 
% figure(7)
% Hcut8=plot_enttropy_sats_CUMULATE(XsigSat_cut8,k,'b',Tvec);

figure(8)
Huterr=plot_est_err_sats_CUMULATE(XsigSat_ut,k,'r',Tvec,ytruth);


figure(9)
Hcut8err=plot_est_err_sats_CUMULATE(XsigSat_cut8,k,'b',Tvec,ytruth);


    pause(0.05)

end

% save('Sat2Dtasking_estimation_sim1','XsigSat_ut','XsigSat_cut8','MeasPairs_ut','MeasPairs_cut8','Hut','Hcut8')

