%% Description of this simulation

% 1. Greedy in time
% 2. Large number of objects
% 3. Make everything parallel as possible



clc
clear
close all

fname='Sat3Dtasking_estimation_greedyTime_';
for simno=1:1000
    if exist(strcat(fname,num2str(simno)),'file')==0
        break
    end
end
simno=simno+1;

load_known_sat_data=false;
saveonoff=1;
saveonoffpic=0;

n=10;
aL=10000;
aU=20000;

%% time step setup
tf=2*24*60*60; % 48 hours
dt=2*60; % 5 mins time step
t0=0;
Tvec=t0:dt:tf;
nt=length(Tvec);
MeasFreq=1; % number of time steps after which measruement update is done

Re=6378.1;
MU=398600.4418;
Tp=@(a)2*pi*sqrt(a^3/MU);

%% satellites

P0=blkdiag(0.1^2,0.1^2,0.1^2,1e-4^2,1e-4^2,1e-4^2);


if load_known_sat_data==true
    load('SatTruth100')
else
    Xsat0=getInitialrv_3D(n,aL,aU);
    
    Nsat=size(Xsat0,1);
    
    
    opt = odeset('reltol',1e-12,'abstol',1e-12);
    ytruth=cell(Nsat,1);
    ytruth_orb=cell(Nsat,1);
    yplottruth=cell(Nsat,1);
    
    parfor i=1:Nsat
        [~,xx]=ode45(@twoBody,Tvec,Xsat0(i,:)',opt);
        ytruth{i}=xx;
        [tt,xx]=ode45(@twoBody,linspace(t0,tf,500),Xsat0(i,:)',opt);
        yplottruth{i}=xx;
        ytruth_orb{i} = XYZ2OE_multiple(ytruth{i});
        i
    end
    
    
end

fn=6;
Nsat=size(Xsat0,1);
Radmodel.Nsat=Nsat;
%% radars (lat,long, altitude) in geod+edic
nsens_out=3;
close all
th_rng=[30*pi/180,330*pi/180];
phi_rng=[-pi/6,pi/6];

[Ang,~]=uniform_sigma_pts([th_rng(1),phi_rng(1)],[th_rng(2),phi_rng(2)],4);


RadPos=zeros(size(Ang,1),3);
MeasCov=zeros(size(Ang,1),nsens_out^2);
k=1;
for i=1:1:size(Ang,1)
    RadPos(k,:)=[Ang(i,1),Ang(i,2),Re];  % [th,phi,Re] Re means on the surface
    R=blkdiag(0.1^2,(0.2*pi/180)^2,(0.2*pi/180)^2);
    %  R=blkdiag((0.1)^2);
    MeasCov(k,:)=reshape(R,1,nsens_out^2);
    k=k+1;
    
end
Nrad=size(RadPos,1);
SensParas=[pi/4*ones(Nrad,1),15000*ones(Nrad,1)];  %misc paras such as [cone_angle,max dist of meas]

Radmodel.SensParas=SensParas;
Radmodel.RadPos=RadPos;
Radmodel.Nrad=Nrad;
Radmodel.hn=nsens_out;
% Radmodel.Q=blkdiag(0.001^2,0.001^2,0.001^2,1e-6^2,1e-6^2,1e-6^2);

Radmodel.G=@(x,Srad)radar_sens_penalty(x,Srad,Radmodel,100);
Radmodel.h=@(x,Srad)radar_sens(x,Srad,Radmodel);

Radmodel.R=@(Srad)reshape(MeasCov(Srad,:),Radmodel.hn,Radmodel.hn);
Radmodel.G=@(x,Srad)radar_sens_penalty(x,Srad,Radmodel,100);

paras.Sathighlight=2;
plot_sat_radar_system2(yplottruth,Radmodel,paras)
xlabel('x')
ylabel('y')
zlabel('z')



%% checking if all the orbits are observable
if load_known_sat_data==false
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
    ndel=find(Satobserve<50);
    if length(ndel)>0
        keyboard
    end
    %     Xsat0(ndel,:)=[];
    %     ytruth(ndel)=[];
    %     yplottruth(ndel)=[];
    %     ytruth_orb(ndel)=[];
    %     Nsat=size(Xsat0,1);
    %
    %      save('SatTruth100','Xsat0','ytruth','yplottruth','ytruth_orb','tf','dt','t0','Tvec','nt','MeasFreq')
    
end






%% Methods that are being used
% UT and CUT8
%generate initial sigma points for each sat
method='ut';
[X,w]=uniform_sigma_pts(zeros(fn,1),eye(fn));
SatState=cell(Nsat,1);
H=zeros(1,Nsat);

for i=1:Nsat
    SatState{i}.mu=zeros(nt,fn);
    SatState{i}.P=zeros(nt,fn^2);
    SatState{i}.Xpts=cell(nt);
    SatState{i}.Wpts=cell(nt);
    SatState{i}.Ymeas_t=zeros(nt,1);
    
    pprr=mvnrnd(Xsat0(i,:),P0);
    SatState{i}.mu(1,:)=pprr;
    SatState{i}.P(1,:)=reshape(P0,1,36);
    H(1,i)=log(det(P0));
    
    SatState{i}.Ymeas_t(1)=1;
    SatState{i}.Q=P0;
    SatState{i}.ytruth=ytruth{i};
end

clear ytruth;



%

%% Propagate, task and filter
MeasPairs=cell(nt,1);  % cell of time steps, the pairing is stored

NTtask=10;
for k=2:1:nt
    % Assuming that everything is done at k-1
    k
    
    
    % Get measuremetn pairs over time from k to k+NTtask
    if rem(k-2,NTtask)==0
        disp(strcat('tasking at ',num2str(k)))
        tic
        if k+NTtask>nt
            MeasPairs=GetMeasPairs_satprob(MeasPairs,SatState,Radmodel,k,nt,Tvec,'ut','GreedyTimeGreedySensors');
        else
            MeasPairs=GetMeasPairs_satprob(MeasPairs,SatState,Radmodel,k,k+NTtask,Tvec,'ut','GreedyTimeGreedySensors');
        end
        toc
    end
    
    
    % Propagation from k-1 to k
    paras.ProcessNoise=true;
    SatState=Propagation_SatState(SatState,k-1,k,Tvec,'ut',paras);
    
    % Measurement update only for the time step k
    paras.TaskAll=false;
    paras.pseudoupdate=false;
    SatState=MeasUpdate_SatState(SatState,MeasPairs,Radmodel,k,Tvec,'ut',paras);
    
    figure(4)
    Hut=plot_enttropy_sats_CUMULATE(SatState,k,'r',Tvec,{'plotit'});
    figure(5)
    Huterr=plot_est_err_sats_CUMULATE_vel(SatState,k,'r',Tvec,ytruth,{'plotit'});
    
    
    %     figure(6)
    %     plot3(ytruth{1}(1:k,1),ytruth{1}(1:k,2),ytruth{1}(1:k,3),'r',ytruth{1}(k,1),ytruth{1}(k,2),ytruth{1}(k,3),'ro','MarkerSize',6)
    %     hold on
    %     plot3(XsigSat_ut{1,1}(1:k,1),XsigSat_ut{1,1}(1:k,2),XsigSat_ut{1,1}(1:k,3),'b',XsigSat_ut{1,1}(k,1),XsigSat_ut{1,1}(k,2),XsigSat_ut{1,1}(k,3),'bo','MarkerSize',6)
    %
    pause(0.05)
    %     keyboard
end


if saveonoff==1
    codeused = fileread('EgSensorTasking_Eg1.m');
    save(strcat(fname,num2str(simno)))
end


%% plotting
close all
for k=507:3:nt
    % Assuming that everything is done at k-1
    k
    
    figure(1)
    clf
    plot_sat_radar_system(XsigSat_ut,k,Radmodel,MeasPairs_ut,yplottruth)
    view([-11,22])
    pause(0.1)
    saveas(gca,strcat(fname,'/sys_',num2str(k)),'png')
    saveas(gca,strcat(fname,'/sys_',num2str(k)),'fig')
    
    %     plot_prop_paper
    
    figure(2)
    clf
    Hut=plot_enttropy_sats_CUMULATE(XsigSat_ut,k,'r',Tvec,{'plotit'});
    plot_prop_paper
    pause(0.1)
    saveas(gca,strcat(fname,'/entr_',num2str(k)),'png')
    saveas(gca,strcat(fname,'/entr_',num2str(k)),'fig')
    
    figure(3)
    clf
    Huterr=plot_est_err_sats_CUMULATE_vel(XsigSat_ut,k,'r',Tvec,ytruth,{'plotit'});
    plot_prop_paper
    pause(0.1)
    saveas(gca,strcat(fname,'/err_',num2str(k)),'png')
    saveas(gca,strcat(fname,'/err_',num2str(k)),'fig')
    
    %     figure(6)
    %     plot3(ytruth{1}(1:k,1),ytruth{1}(1:k,2),ytruth{1}(1:k,3),'r',ytruth{1}(k,1),ytruth{1}(k,2),ytruth{1}(k,3),'ro','MarkerSize',6)
    %     hold on
    %     plot3(XsigSat_ut{1,1}(1:k,1),XsigSat_ut{1,1}(1:k,2),XsigSat_ut{1,1}(1:k,3),'b',XsigSat_ut{1,1}(k,1),XsigSat_ut{1,1}(k,2),XsigSat_ut{1,1}(k,3),'bo','MarkerSize',6)
    %
    %     pause(0.1)
    %     keyboard
end


%% Analytics
A=cell(1,2);
jj=1;
for k=1:1:nt
    A{1,k+2}=Tvec(k)/3600;
end
for i=1:1:Nsat
    A{2*i,1}=i;
    pp=RVtoCOEs(Xsat0(i,1:3),Xsat0(i,4:6));
    A{2*i,2}=Tp(pp(1))/3600;
    for k=1:1:nt
        P=XsigSat_ut{i,2}(k,:);
        n=sqrt(length(P));
        A{2*i+1,k+2}=sqrt(max(eig(reshape(P,n,n))));
        
        A{2*i,k+2}=-1;
        for nrad=1:1:Nrad
            [yy,hh]=Radmodel.G(ytruth{i}(k,:),nrad);
            if isnan(hh)==0
                A{2*i,k+2}=0;
            end
        end
        if isempty(MeasPairs_ut{k})==0
            if length(find(MeasPairs_ut{k}(:,1)==i))>0 % it has been paired
                A{2*i,k+2}=1;
            end
        end
        
    end
end
cell2csv('TaskingDebug.csv',A,',')

