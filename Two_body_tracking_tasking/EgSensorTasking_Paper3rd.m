%eg 3: multiple sensor and multiple moving target with coverage problem
% Initiating the example
clc
close all
clear all
t0=0;
dt=1;
tf=100;
TC=TargSens(t0,dt,tf);

 TC.SimProps.Xboundary=[0 100];
 TC.SimProps.Yboundary=[0 100];
 
disp('done initializing')
TC=TC.Set_Sigmapts('UT');
disp('done Sig pt setup')
% Adding the sensors

% add the tracking sensors
% TC=TC.Add_Sensor(obj,SensorType,SensorModel,SensorConstraints,R,sk,alphak,rmaxk,dirn,FOVpenaltytype)
% TC=TC.Add_Sensor('Stationary','Range+Bearing','1Sensor->1Target',diag([1^2,(0.5*pi/180)^2]),[10,10],pi/8,30,pi/6,'Simple');
% TC=TC.Add_Sensor('Stationary','Range+Bearing','1Sensor->1Target',diag([2^2,(1*pi/180)^2]),[90,90],pi/8,60,-pi/6,'Simple');

%add the tracking/coverage sensors
% TC=TC.Add_MultiModalSensor(obj,nmodes,SensorType,SensorModel,SensorConstraints,R,sk,alphak,rmaxk,dirn,FOVpenaltytype)

 
SENSCONSTRAINTS.Seq={};
SENSCONSTRAINTS.deq={};
SENSCONSTRAINTS.Sineq={};
SENSCONSTRAINTS.dineq={};         

% #############################   SENSOR 1  ###############################



sens.SensorsID=1;
sens.SensorDynTypes_allk=repmat({'Stationary'},TC.SimProps.Time.nsteps,1);
sens.SensorNmodes_allk=repmat(2,TC.SimProps.Time.nsteps,1);
sens.R_allk= repmat({diag([0.2^2,(0.1*pi/180)^2])},TC.SimProps.Time.nsteps,1); 
sens.SensorModel_allk=repmat({'Range+Bearing'},TC.SimProps.Time.nsteps,1); 
sens.SensorFOVpenalty_allk=repmat({'Simple'},TC.SimProps.Time.nsteps,1); 
sens.SensorFOVtype_allk=repmat({'1Sensor->1Target'},TC.SimProps.Time.nsteps,1); 
sens.SensorDynMotionModel_allk=repmat({'None'},TC.SimProps.Time.nsteps,1); 
sens.SensorPos_allk=repmat([30,30],TC.SimProps.Time.nsteps,1); 
sens.SensorRmax_allk=repmat(70,TC.SimProps.Time.nsteps,1); 
sens.SensorDir_allk=repmat(0,TC.SimProps.Time.nsteps,1); 
sens.SensorAlpha_allk=repmat(pi/10,TC.SimProps.Time.nsteps,1); 
sens.SensorModeConstraints_allk=repmat({'1ModeOperation'},TC.SimProps.Time.nsteps,1); 
sens.SensorOnOff_allk=repmat({'On'},TC.SimProps.Time.nsteps,1); 

TC=TC.Add_Sensor(sens);

sens.SensorsID=2;
sens.SensorDynTypes_allk=repmat({'Stationary'},TC.SimProps.Time.nsteps,1);
sens.SensorNmodes_allk=repmat(2,TC.SimProps.Time.nsteps,1);
sens.R_allk= repmat({diag([1^2,(1*pi/180)^2])},TC.SimProps.Time.nsteps,1); 
sens.SensorModel_allk=repmat({'Range+Bearing'},TC.SimProps.Time.nsteps,1); 
sens.SensorFOVpenalty_allk=repmat({'Simple'},TC.SimProps.Time.nsteps,1); 
sens.SensorFOVtype_allk=repmat({'1Sensor->1Target'},TC.SimProps.Time.nsteps,1); 
sens.SensorDynMotionModel_allk=repmat({'None'},TC.SimProps.Time.nsteps,1); 
sens.SensorPos_allk=repmat([30,30],TC.SimProps.Time.nsteps,1); 
sens.SensorRmax_allk=repmat(50,TC.SimProps.Time.nsteps,1); 
sens.SensorDir_allk=repmat(0,TC.SimProps.Time.nsteps,1); 
sens.SensorAlpha_allk=repmat(pi,TC.SimProps.Time.nsteps,1); 
sens.SensorModeConstraints_allk=repmat({'1ModeOperation'},TC.SimProps.Time.nsteps,1); 
sens.SensorOnOff_allk=repmat({'On'},TC.SimProps.Time.nsteps,1); 

TC=TC.Add_Sensor(sens);

SENSCONSTRAINTS.Seq=vertcat(SENSCONSTRAINTS.Seq,{[1,2]});
SENSCONSTRAINTS.deq=vertcat(SENSCONSTRAINTS.deq,{1});

% #############################   SENSOR 2  ###############################

sens.SensorsID=3;
sens.SensorDynTypes_allk=repmat({'Stationary'},TC.SimProps.Time.nsteps,1);
sens.SensorNmodes_allk=repmat(2,TC.SimProps.Time.nsteps,1);
sens.R_allk= repmat({diag([0.2^2,(0.1*pi/180)^2])},TC.SimProps.Time.nsteps,1); 
sens.SensorModel_allk=repmat({'Range+Bearing'},TC.SimProps.Time.nsteps,1); 
sens.SensorFOVpenalty_allk=repmat({'Simple'},TC.SimProps.Time.nsteps,1); 
sens.SensorFOVtype_allk=repmat({'1Sensor->1Target'},TC.SimProps.Time.nsteps,1); 
sens.SensorDynMotionModel_allk=repmat({'None'},TC.SimProps.Time.nsteps,1); 
sens.SensorPos_allk=repmat([70,70],TC.SimProps.Time.nsteps,1); 
sens.SensorRmax_allk=repmat(70,TC.SimProps.Time.nsteps,1); 
sens.SensorDir_allk=repmat(0,TC.SimProps.Time.nsteps,1); 
sens.SensorAlpha_allk=repmat(pi/10,TC.SimProps.Time.nsteps,1); 
sens.SensorModeConstraints_allk=repmat({'1ModeOperation'},TC.SimProps.Time.nsteps,1); 
sens.SensorOnOff_allk=repmat({'On'},TC.SimProps.Time.nsteps,1); 

TC=TC.Add_Sensor(sens);

sens.SensorsID=4;
sens.SensorDynTypes_allk=repmat({'Stationary'},TC.SimProps.Time.nsteps,1);
sens.SensorNmodes_allk=repmat(2,TC.SimProps.Time.nsteps,1);
sens.R_allk= repmat({diag([1^2,(1*pi/180)^2])},TC.SimProps.Time.nsteps,1); 
sens.SensorModel_allk=repmat({'Range+Bearing'},TC.SimProps.Time.nsteps,1); 
sens.SensorFOVpenalty_allk=repmat({'Simple'},TC.SimProps.Time.nsteps,1); 
sens.SensorFOVtype_allk=repmat({'1Sensor->1Target'},TC.SimProps.Time.nsteps,1); 
sens.SensorDynMotionModel_allk=repmat({'None'},TC.SimProps.Time.nsteps,1); 
sens.SensorPos_allk=repmat([70,70],TC.SimProps.Time.nsteps,1); 
sens.SensorRmax_allk=repmat(50,TC.SimProps.Time.nsteps,1); 
sens.SensorDir_allk=repmat(0,TC.SimProps.Time.nsteps,1); 
sens.SensorAlpha_allk=repmat(pi,TC.SimProps.Time.nsteps,1); 
sens.SensorModeConstraints_allk=repmat({'1ModeOperation'},TC.SimProps.Time.nsteps,1); 
sens.SensorOnOff_allk=repmat({'On'},TC.SimProps.Time.nsteps,1); 

TC=TC.Add_Sensor(sens);

SENSCONSTRAINTS.Seq=vertcat(SENSCONSTRAINTS.Seq,{[3,4]});
SENSCONSTRAINTS.deq=vertcat(SENSCONSTRAINTS.deq,{1});

% #############################   SENSOR 3  ###############################


sens.SensorsID=5;
sens.SensorDynTypes_allk=repmat({'Stationary'},TC.SimProps.Time.nsteps,1);
sens.SensorNmodes_allk=repmat(2,TC.SimProps.Time.nsteps,1);
sens.R_allk= repmat({diag([0.2^2,(0.1*pi/180)^2])},TC.SimProps.Time.nsteps,1); 
sens.SensorModel_allk=repmat({'Range+Bearing'},TC.SimProps.Time.nsteps,1); 
sens.SensorFOVpenalty_allk=repmat({'Simple'},TC.SimProps.Time.nsteps,1); 
sens.SensorFOVtype_allk=repmat({'1Sensor->1Target'},TC.SimProps.Time.nsteps,1); 
sens.SensorDynMotionModel_allk=repmat({'None'},TC.SimProps.Time.nsteps,1); 
sens.SensorPos_allk=repmat([70,30],TC.SimProps.Time.nsteps,1); 
sens.SensorRmax_allk=repmat(70,TC.SimProps.Time.nsteps,1); 
sens.SensorDir_allk=repmat(0,TC.SimProps.Time.nsteps,1); 
sens.SensorAlpha_allk=repmat(pi/10,TC.SimProps.Time.nsteps,1); 
sens.SensorModeConstraints_allk=repmat({'1ModeOperation'},TC.SimProps.Time.nsteps,1); 
sens.SensorOnOff_allk=repmat({'On'},TC.SimProps.Time.nsteps,1); 

TC=TC.Add_Sensor(sens);

sens.SensorsID=6;
sens.SensorDynTypes_allk=repmat({'Stationary'},TC.SimProps.Time.nsteps,1);
sens.SensorNmodes_allk=repmat(2,TC.SimProps.Time.nsteps,1);
sens.R_allk= repmat({diag([1^2,(1*pi/180)^2])},TC.SimProps.Time.nsteps,1); 
sens.SensorModel_allk=repmat({'Range+Bearing'},TC.SimProps.Time.nsteps,1); 
sens.SensorFOVpenalty_allk=repmat({'Simple'},TC.SimProps.Time.nsteps,1); 
sens.SensorFOVtype_allk=repmat({'1Sensor->1Target'},TC.SimProps.Time.nsteps,1); 
sens.SensorDynMotionModel_allk=repmat({'None'},TC.SimProps.Time.nsteps,1); 
sens.SensorPos_allk=repmat([70,30],TC.SimProps.Time.nsteps,1); 
sens.SensorRmax_allk=repmat(50,TC.SimProps.Time.nsteps,1); 
sens.SensorDir_allk=repmat(0,TC.SimProps.Time.nsteps,1); 
sens.SensorAlpha_allk=repmat(pi,TC.SimProps.Time.nsteps,1); 
sens.SensorModeConstraints_allk=repmat({'1ModeOperation'},TC.SimProps.Time.nsteps,1); 
sens.SensorOnOff_allk=repmat({'On'},TC.SimProps.Time.nsteps,1); 

TC=TC.Add_Sensor(sens);

SENSCONSTRAINTS.Seq=vertcat(SENSCONSTRAINTS.Seq,{[5,6]});
SENSCONSTRAINTS.deq=vertcat(SENSCONSTRAINTS.deq,{1});

% #############################   SENSOR 4  ###############################


sens.SensorsID=7;
sens.SensorDynTypes_allk=repmat({'Stationary'},TC.SimProps.Time.nsteps,1);
sens.SensorNmodes_allk=repmat(2,TC.SimProps.Time.nsteps,1);
sens.R_allk= repmat({diag([0.2^2,(0.1*pi/180)^2])},TC.SimProps.Time.nsteps,1); 
sens.SensorModel_allk=repmat({'Range+Bearing'},TC.SimProps.Time.nsteps,1); 
sens.SensorFOVpenalty_allk=repmat({'Simple'},TC.SimProps.Time.nsteps,1); 
sens.SensorFOVtype_allk=repmat({'1Sensor->1Target'},TC.SimProps.Time.nsteps,1); 
sens.SensorDynMotionModel_allk=repmat({'None'},TC.SimProps.Time.nsteps,1); 
sens.SensorPos_allk=repmat([30,70],TC.SimProps.Time.nsteps,1); 
sens.SensorRmax_allk=repmat(70,TC.SimProps.Time.nsteps,1); 
sens.SensorDir_allk=repmat(0,TC.SimProps.Time.nsteps,1); 
sens.SensorAlpha_allk=repmat(pi/10,TC.SimProps.Time.nsteps,1); 
sens.SensorModeConstraints_allk=repmat({'1ModeOperation'},TC.SimProps.Time.nsteps,1); 
sens.SensorOnOff_allk=repmat({'On'},TC.SimProps.Time.nsteps,1); 

TC=TC.Add_Sensor(sens);

sens.SensorsID=8;
sens.SensorDynTypes_allk=repmat({'Stationary'},TC.SimProps.Time.nsteps,1);
sens.SensorNmodes_allk=repmat(2,TC.SimProps.Time.nsteps,1);
sens.R_allk= repmat({diag([1^2,(1*pi/180)^2])},TC.SimProps.Time.nsteps,1); 
sens.SensorModel_allk=repmat({'Range+Bearing'},TC.SimProps.Time.nsteps,1); 
sens.SensorFOVpenalty_allk=repmat({'Simple'},TC.SimProps.Time.nsteps,1); 
sens.SensorFOVtype_allk=repmat({'1Sensor->1Target'},TC.SimProps.Time.nsteps,1); 
sens.SensorDynMotionModel_allk=repmat({'None'},TC.SimProps.Time.nsteps,1); 
sens.SensorPos_allk=repmat([30,70],TC.SimProps.Time.nsteps,1); 
sens.SensorRmax_allk=repmat(50,TC.SimProps.Time.nsteps,1); 
sens.SensorDir_allk=repmat(0,TC.SimProps.Time.nsteps,1); 
sens.SensorAlpha_allk=repmat(pi,TC.SimProps.Time.nsteps,1); 
sens.SensorModeConstraints_allk=repmat({'1ModeOperation'},TC.SimProps.Time.nsteps,1); 
sens.SensorOnOff_allk=repmat({'On'},TC.SimProps.Time.nsteps,1); 

TC=TC.Add_Sensor(sens);

SENSCONSTRAINTS.Seq=vertcat(SENSCONSTRAINTS.Seq,{[7,8]});
SENSCONSTRAINTS.deq=vertcat(SENSCONSTRAINTS.deq,{1});

TC=TC.UpdateSensorConstraints(SENSCONSTRAINTS);

disp('done Sensor Placement')


%%  Adding targets

% Adding real targets
% TC=TC.Add_Target(TargType,DynType,xk,Pk,truth,Qk);

%============= target 1 ===================================================

xktruth=zeros(TC.SimProps.Time.nsteps,5);
xktruth(1,:)=[1,1,3,3.5,0];
for t=2:1:TC.SimProps.Time.nsteps
    xktruth(t,:)=KIRB_CT_eg_dyn_disc(xktruth(t-1,:)',dt);
    xktruth(t,5)=sin(t/3);  % give some turn rate to get some interesting trajectory
end

Pk=zeros(TC.SimProps.Time.nsteps,25);
Pk(1,:)=reshape(diag([3^2,3^2,0.9^2,0.9^2,(1*pi/180)^2]),1,25);
xk=zeros(TC.SimProps.Time.nsteps,5);
xk(1,:)=xktruth(1,:);


targ.TargetID=1;
targ.MeanState_allk=xk;
targ.CovState_allk=Pk;
targ.VisibilityState_allk=repmat({'Tracked'},TC.SimProps.Time.nsteps,1);  
targ.TruePosState_allk=xktruth;
targ.TargetDynTypes_allk=repmat({'Move'},TC.SimProps.Time.nsteps,1);  
targ.Q_allk=repmat({diag([2^2,2^2,0.01^2,0.01^2,(0.1*pi/180)^2])},TC.SimProps.Time.nsteps,1); 
targ.TargetModel_allk=repmat({'CT'},TC.SimProps.Time.nsteps,1);
targ.TargetModelDim_allk=repmat(5,TC.SimProps.Time.nsteps,1);
            
TC=TC.Add_Target(targ);

%============= target 2 ===================================================

xktruth=zeros(TC.SimProps.Time.nsteps,5);
xktruth(1,:)=[15,90,-3,-3.5,0];
for t=2:1:TC.SimProps.Time.nsteps
    xktruth(t,:)=KIRB_CT_eg_dyn_disc(xktruth(t-1,:)',dt);
    xktruth(t,5)=sin(t/2);  % give some turn rate to get some interesting trajectory
end
Pk=zeros(TC.SimProps.Time.nsteps,25);
Pk(1,:)=reshape(diag([3^2,3^2,0.9^2,0.9^2,(1*pi/180)^2]),1,25);
xk=zeros(TC.SimProps.Time.nsteps,5);
xk(1,:)=xktruth(1,:);


targ.TargetID=2;
targ.MeanState_allk=xk;
targ.CovState_allk=Pk;
targ.VisibilityState_allk=repmat({'Tracked'},TC.SimProps.Time.nsteps,1);  
targ.TruePosState_allk=xktruth;
targ.TargetDynTypes_allk=repmat({'Move'},TC.SimProps.Time.nsteps,1);  
targ.Q_allk=repmat({diag([2^2,2^2,0.01^2,0.01^2,(0.1*pi/180)^2])},TC.SimProps.Time.nsteps,1); 
targ.TargetModel_allk=repmat({'CT'},TC.SimProps.Time.nsteps,1);
targ.TargetModelDim_allk=repmat(5,TC.SimProps.Time.nsteps,1);
      
TC=TC.Add_Target(targ);

%============= target 3 ===================================================

xktruth=zeros(TC.SimProps.Time.nsteps,4);
xktruth(1,:)=[80,90,-3,-0];
for t=2:1:TC.SimProps.Time.nsteps
    xktruth(t,:)=KIRB_UM_eg_dyn_disc(xktruth(t-1,:)',dt);
    if rem(t,10)==0
        th=pi/4;
        xktruth(t,3:4)=[cos(th),-sin(th);sin(th),cos(th)]*xktruth(t,3:4)';  % give some turn rate to get some interesting trajectory
    end
end
Pk=zeros(TC.SimProps.Time.nsteps,16);
Pk(1,:)=reshape(diag([3^2,3^2,0.9^2,0.9^2]),1,25);
xk=zeros(TC.SimProps.Time.nsteps,4);
xk(1,:)=xktruth(1,:);



targ.TargetID=3;
targ.MeanState_allk=xk;
targ.CovState_allk=Pk;
targ.VisibilityState_allk=repmat({'Tracked'},TC.SimProps.Time.nsteps,1);  
targ.TruePosState_allk=xktruth;
targ.TargetDynTypes_allk=repmat({'Move'},TC.SimProps.Time.nsteps,1);  
targ.Q_allk=repmat({diag([2^2,2^2,0.01^2,0.01^2])},TC.SimProps.Time.nsteps,1); 
targ.TargetModel_allk=repmat({'UM'},TC.SimProps.Time.nsteps,1);
targ.TargetModelDim_allk=repmat(4,TC.SimProps.Time.nsteps,1);

TC=TC.Add_Target(targ);

% %%  adding virtual objects
% %% ============= targets v1 ===================================================
[Xvirtgrid,Yvirtgrid]=meshgrid(TC.SimProps.Xboundary(1):50:TC.SimProps.Xboundary(2),TC.SimProps.Yboundary(1):50:TC.SimProps.Yboundary(2));
kk=1;
for i=1:1:size(Xvirtgrid,1)
    for j=1:1:size(Xvirtgrid,2)
        
        xktruth=zeros(TC.SimProps.Time.nsteps,2);
        xktruth(1,:)=[Xvirtgrid(i,j),Yvirtgrid(i,j)];
        for t=2:1:TC.SimProps.Time.nsteps
            xktruth(t,:)=xktruth(t-1,:);
        end
        
        Pk=zeros(TC.SimProps.Time.nsteps,4);
        Pk(1,:)=reshape(diag([3^2,3^2]),1,4);
        xk=zeros(TC.SimProps.Time.nsteps,2);
        xk(1,:)=xktruth(1,:);
        

        targ.TargetID=3+kk;
        kk=kk+1;
        targ.MeanState_allk=xk;
        targ.CovState_allk=Pk;
        targ.VisibilityState_allk=repmat({'Tracked'},TC.SimProps.Time.nsteps,1);
        targ.TruePosState_allk=xktruth;
        targ.TargetDynTypes_allk=repmat({'VirtualStationary'},TC.SimProps.Time.nsteps,1);
        targ.Q_allk=repmat({diag([2^2,2^2])},TC.SimProps.Time.nsteps,1);
        targ.TargetModel_allk=repmat({'NoDyn'},TC.SimProps.Time.nsteps,1);
        targ.TargetModelDim_allk=repmat(2,TC.SimProps.Time.nsteps,1);
        
        TC=TC.Add_Target(targ);
    end
end



disp('done target Setup')
%% Initial Debug plot

TC.plot_truth_scenario()
disp('done initial plot')

TC.PrintDebug()

%% Debug tasking
TC=TC.OptiSens_TargPairs('MIUB');
TC.plot_Dynamic_scenario()

%%  Propagating the targets

TC=TC.SimulationModeOption(0);

for k=2:1:TC.nt
TC=TC.DynProp_Sensor(TC.tk,1:TC.Ns); % from tk to tk+1... without change in tk
TC=TC.DynProp_Target(1:TC.Nt);       % from tk to tk+1... without change in tk


TC=TC.UpdateTimek(); %update to tk+1 i.e. tk becomes tk+1
TC=TC.OptiSens_TargPairs('MIUB'); % do sensor tasking at the new tk (i.e. tk+1)
TC.plot_Dynamic_scenario();

TC.FIM
TC.MU
keyboard

TC=TC.MeasUpdate_Target(1:TC.Nt); % Do the measurement update of all the targets


% pause

disp('              ')
disp(num2str(k))
disp('              ')
end
