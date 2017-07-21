%eg 3: multiple sensor and multiple moving target with coverage problem
%% Initiating the example
clc
close all
clear all
t0=0;
dt=1;
tf=100;
TC=TargSens(t0,dt,tf);

TC.xlim=100;
TC.ylim=100;

disp('done initializing')
TC=TC.Set_Sigmapts('UT');
disp('done Sig pt setup')
%% Adding the sensors

%add the tracking sensors
% TC=TC.Add_Sensor(obj,SensorType,SensorModel,SensorConstraints,R,sk,alphak,rmaxk,dirn,FOVpenaltytype)
% TC=TC.Add_Sensor('Stationary','Range+Bearing','1Sensor->1Target',diag([1^2,(0.5*pi/180)^2]),[10,10],pi/8,30,pi/6,'Simple');
% TC=TC.Add_Sensor('Stationary','Range+Bearing','1Sensor->1Target',diag([2^2,(1*pi/180)^2]),[90,90],pi/8,60,-pi/6,'Simple');


%add the tracking/coverage sensors
% TC=TC.Add_MultiModalSensor(obj,nmodes,SensorType,SensorModel,SensorConstraints,R,sk,alphak,rmaxk,dirn,FOVpenaltytype)
TC=TC.Add_MultiModalSensor(2,{'Stationary','Stationary'},{'Range+Bearing','Range+Bearing'},...
    {'1Sensor->1Target','1Sensor->AllTarget'},{diag([0.2^2,(0.2*pi/180)^2]),diag([2^2,(2*pi/180)^2])}...
    ,[30,30],{pi/10,pi},{70,50},{0,0},'Simple');
TC=TC.Add_MultiModalSensor(2,{'Stationary','Stationary'},{'Range+Bearing','Range+Bearing'},...
    {'1Sensor->1Target','1Sensor->AllTarget'},{diag([0.2^2,(0.2*pi/180)^2]),diag([2^2,(2*pi/180)^2])}...
    ,[70,70],{pi/10,pi},{70,50},{pi/4,0},'Simple');
TC=TC.Add_MultiModalSensor(2,{'Stationary','Stationary'},{'Range+Bearing','Range+Bearing'},...
    {'1Sensor->1Target','1Sensor->AllTarget'},{diag([0.2^2,(0.2*pi/180)^2]),diag([2^2,(2*pi/180)^2])}...
    ,[70,30],{pi/10,pi},{70,50},{pi/4,0},'Simple');
TC=TC.Add_MultiModalSensor(2,{'Stationary','Stationary'},{'Range+Bearing','Range+Bearing'},...
    {'1Sensor->1Target','1Sensor->AllTarget'},{diag([0.2^2,(0.2*pi/180)^2]),diag([2^2,(2*pi/180)^2])}...
    ,[30,70],{pi/10,pi},{70,50},{pi/4,0},'Simple');

disp('done Sensor Placement')
%% Adding targets
% Adding real targets
% TC=TC.Add_Target(TargType,DynType,xk,Pk,truth,Qk);

%============= target 1 ===================================================
xktruth=zeros(TC.nt,5);
xktruth(1,:)=[1,1,3,3.5,0];
for t=2:1:TC.nt
    xktruth(t,:)=KIRB_CT_eg_dyn_disc(xktruth(t-1,:)',dt);
    xktruth(t,5)=sin(t/3);  % give some turn rate to get some interesting trajectory
end
Pk=diag([2^2,2^2,0.9^2,0.9^2,(1*pi/180)^2]);
xk=mvnrnd(xktruth(1,:),Pk);
Qk=diag([0.2^2,0.2^2,0.01^2,0.01^2,(0.1*pi/180)^2]);
TC=TC.Add_Target('Move','CT',xk,Pk,xktruth,Qk);

%============= target 2 ===================================================
xktruth=zeros(TC.nt,5);
xktruth(1,:)=[15,90,-3,-3.5,0];
for t=2:1:TC.nt
    xktruth(t,:)=KIRB_CT_eg_dyn_disc(xktruth(t-1,:)',dt);
    xktruth(t,5)=sin(t/2);  % give some turn rate to get some interesting trajectory
end
Pk=diag([2^2,2^2,0.9^2,0.9^2,(1*pi/180)^2]);
xk=mvnrnd(xktruth(1,:),Pk);
Qk=diag([0.2^2,0.2^2,0.01^2,0.01^2,(0.1*pi/180)^2]);
TC=TC.Add_Target('Move','CT',xk,Pk,xktruth,Qk);

%============= target 3===================================================
xktruth=zeros(TC.nt,4);
xktruth(1,:)=[80,90,-3,-0];
for t=2:1:TC.nt
    xktruth(t,:)=KIRB_UM_eg_dyn_disc(xktruth(t-1,:)',dt);
    if rem(t,10)==0
        th=pi/4;
        xktruth(t,3:4)=[cos(th),-sin(th);sin(th),cos(th)]*xktruth(t,3:4)';  % give some turn rate to get some interesting trajectory
    end
end
Pk=diag([2^2,2^2,0.9^2,0.9^2]);
xk=mvnrnd(xktruth(1,:),Pk);
Qk=diag([0.2^2,0.2^2,0.01^2,0.01^2]);
TC=TC.Add_Target('Move','UM',xk,Pk,xktruth,Qk);

%adding virtual objects
%============= targets v1===================================================
[Xvirtgrid,Yvirtgrid]=meshgrid(0:15:TC.xlim,0:15:TC.ylim);
for i=1:1:size(Xvirtgrid,1)
    for j=1:1:size(Xvirtgrid,2)
        xktruth=zeros(TC.nt,2);
        xktruth(1,:)=[Xvirtgrid(i,j),Yvirtgrid(i,j)];
        for t=2:1:TC.nt
            xktruth(t,:)=xktruth(t-1,:);
        end
        Pk=diag([2^2,2^2]);
        xk=mvnrnd(xktruth(1,:),Pk);
        Qk=diag([0.2^2,0.2^2]);
        TC=TC.Add_Target('Virtual','NoDyn',xk,Pk,xktruth,Qk);
    end
end


disp('done target Setup')
%% Initial Debug plot

TC.plot_truth_scenario()
disp('done initial plot')


%% Debug tasking
TC.OptiSens_TargPairs('MIUB')


%%







