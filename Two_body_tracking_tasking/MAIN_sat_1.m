% Satellite examples,
% 1.First Construct a scenario 15 satellites and 5,000 satellites
% 2. Do linear all exclusive mutual information --> over
% time,sensors,targets ----->  Take in all the covariances and solve for
% the assignment variables
% 2. Exhaustive Search: Greedy time, Greedy Sensor, Greedy Target ----->  Take in all the covariances and solve for
% the assignment variables

%% Simulation options
close all

redoSATTRAJ=0;
%% Constants

Nsat=15;
Nrad=11;
dt=5*60; % in seconds
Tf=(24)*60*60; % in seconds
Tvec=0:dt:Tf;
plotTvec=0:5*60:Tf;
Ntimesteps=length(Tvec);
Re=6378.1;

SensTaskHorizon=10;

%% Simulation Structures

Constants.Nsat=Nsat;
Constants.Nrad=Nrad;
Constants.dt=dt;
Constants.Tf=Tf;
Constants.Tvec=Tvec;
Constants.plotTvec=plotTvec;
Constants.Ntimesteps=Ntimesteps;
Constants.Re=Re;
Constants.SensTaskHorizon=SensTaskHorizon;

Radars=cell(Nrad,1);
Satellites=cell(Nsat,1);


MeasPairs=cell(Constants.Ntimesteps ,1);
for i=1:1:Constants.Ntimesteps
    MeasPairs{i}=zeros(Constants.Nsat,Constants.Nrad);
end
% MeasPairs{1}(i,j)
%               j=1,rad1  j=2,rad2   j=3,rad3  j=4,rad4   ...
% i=1 sat1
% i=2 sat2
% i=3 sat3
% .
% .
% .



%% radars (lat,long, altitude) in geod+edic
nsens_out=2;
% (th,phi,ConeAngle,MaxRange) --> (perp to equaltor, along equator)
Ang=[-100*pi/180    0  pi/3 2500
    -80*pi/180    20*pi/180 pi/3 2500
    -80*pi/180    -20*pi/180 pi/3 2500
    
    115*pi/180    90*pi/180 pi/3 2500
    110*pi/180    0*pi/180 pi/3 2500
    80*pi/180    0*pi/180 pi/3 2500
    75*pi/180    90*pi/180 pi/3 2500
    
    0.1437   -0.0168 pi/3 1500
    -0.1357   -0.0114 pi/3 1500
    
    0.1663   -3.0993 pi/3 3000
    -0.1549   -3.0988 pi/3 3000];

RadPos=zeros(size(Ang,1),3);
k=1;

for i=1:1:Constants.Nrad
    Radars{i}.PolarPositions=[];
    Radars{i}.R =[];
    Radars{i}.hn=nsens_out;
    Radars{i}.G=@(x,Srad)radar_sens_penalty(x,Srad,Radmodel,100);
    Radars{i}.h=@(x,Srad)radar_sens(x,Srad,Radmodel);
    
    %
    Radars{i}.PolarPositions=vertcat(Radars{i}.PolarPositions,[Ang(i,1),Ang(i,2),Constants.Re]);
    Radars{i}.ConeAngle=Ang(i,3);
    Radars{i}.MaxRange=Ang(i,4);
    
    Radars{i}.R=vertcat(Radars{i}.R, reshape(R,1,Radars{i}.hn^2) );
end




%% Set Sattelite Properties
% figure
% hold on
P0=blkdiag(0.01,0.01,0.01,1e-8,1e-8,1e-8);

for i=1:1:Constants.Nsat
    Satellites{i}.HighlightPlotTraj=0 ;
    Satellites{i}.Q=zeros(6,6);
end

if redoSATTRAJ==1
    Xsat0=getInitialrv_3D_3(Constants.Nsat);
    opt = odeset('reltol',1e-12,'abstol',1e-12);
    ytruth=cell(Nsat,1);
    ytruth_orb=cell(Nsat,1);
    yplottruth=cell(Nsat,1);
    
    
    parfor i=1:Nsat
        
        [X,w]=conjugate_dir_gausspts_till_8moment(Xsat0(i,:)',P0);
        [a,b]=max(sqrt(sum((X(:,4:6)-repmat(Xsat0(i,4:6),length(w),1)).^2,2)));
        pprr=X(b,:);
        
        [~,xx]=ode45(@twoBody,Tvec,pprr',opt);
        ytruth{i}=xx;
        Satellites{i}.ytruth=xx;
        
        [~,xx]=ode45(@twoBody,plotTvec,pprr',opt);
        yplottruth{i}=xx;
        Satellites{i}.yplottruth=xx;
        
        %     ytruth_orb{i} = XYZ2OE_multiple(ytruth{i});
        i
        % plot3(xx(:,1),xx(:,2),xx(:,3))
        % keyboard
    end
    save('SavedData/SATELLIIETASKING_MAIN_sat_1__INITIALSATtrajs.mat','Xsat0','ytruth','yplottruth')
else
    M=load('SavedData/SATELLIIETASKING_MAIN_sat_1__INITIALSATtrajs.mat');
    Xsat0=M.Xsat0;
    ytruth=M.ytruth;
    yplottruth=M.yplottruth;
    
    for i=1:Constants.Nsat
        Satellites{i}.ytruth=ytruth{i};
        Satellites{i}.yplottruth=yplottruth{i};
    end
    
end




disp('print done sat prop')




%% Plot trajectories to verify

plot_sat_radar_system2(Satellites,Radars,Constants)


%% Running Simmulation
NextSensTaskTimeStep=2;

for k=2:1:nt
    disp(strcat('at time step : ',num2str(k)))
    
    tic
    % propagate
    Satellites=Propagate_sattask(Satellites,Constants,k-1,k,'ut');
    
    %Generate the MeasPairs
    if k==NextSensTaskTimeStep
        MeasPairs=SensorTask_FIM(MeasPairs,Satellites,Radars,Constants,k,k+Constants.SensTaskHorizon,'ut');
        NextSensTaskTimeStep=k+Constants.SensTaskHorizon;
    end
    
    % Measurement update only for the MeasPairs{k}
    Satellites=MeasUpdate_sattask(MeasPairs,Satellites,Radars,Constants,k,'ut');
    
    toc
    
    
    pause(0.05)
end
