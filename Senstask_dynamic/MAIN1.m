
% MAIN function to run dynamic sensing
% First build a grid: UAV can jump to a limited number of grid points
% around its current position
% We have to do greedy time, greedy sensor and greedy target.
% Compute the joint MI over time

%% Simulation options
close all
clear all
clc

ss=datestr(now());
filename=['SavedData/UAV_test',ss,'.mat'];

%%

Targets.Ntargs=10;
Sensors.Nsens=2;

time.dt=5; % in seconds
time.Tf=250; % in seconds
time.Tvec=0:time.dt:time.Tf;
time.Ntimesteps=length(time.Tvec);

time.SensTaskHorizon=3; % +1 steps in all

Targets.xf=cell(Targets.Ntargs,1);
Targets.Pf=cell(Targets.Ntargs,1);

Sensors.xfgridpos=cell(Sensors.Nsens,1);
Sensors.xf=cell(Sensors.Nsens,1);



Xlim=[0,500];
dx=50;
Grid=SetupUAV_Grid(Xlim,dx);

%% targets and UAV sensors setup
close all
Targets.f=cell(Targets.Ntargs,1);
Targets.fn=zeros(Targets.Ntargs,1);
Targets.truth=cell(Targets.Ntargs,1);
for i=1:5
    Targets.f{i}=@Stat_location;
    Targets.fn(i)=2;
end
for i=6:10
    Targets.f{i}=@UM_motion;
    Targets.fn(i)=4;
end
% for i=9:10
%     Targets.f{i}=@CT_motion;
%     Targets.fn(i)=5;
% end

% first stationary targets
for i=1:5
    Targets.truth{i}=zeros(time.Ntimesteps,Targets.fn(i));
    p=linspace(Xlim(1),Xlim(2),10);
    a=p(2);
    b=p(end-1);
    x0=[a+(b-a)*rand,a+(b-a)*rand];
    Targets.truth{i}=repmat(x0,time.Ntimesteps,1);
    
    P0 = diag([10^2,10^2]);
    Targets.Q{i}=0.01*diag([1^2,1^2]);
    Targets.xf{i}=zeros(time.Ntimesteps,Targets.fn(i));
    Targets.Pf{i}=zeros(time.Ntimesteps,Targets.fn(i)^2);
    Targets.Pf{i}(1,:) =  reshape(P0,1,Targets.fn(i)^2);
    Targets.xf{i}(1,:) = mvnrnd(Targets.truth{i}(1,:),P0);
    
end

% uniform motion targets
mindist = 400;
for i=6:10
    Targets.truth{i}=zeros(time.Ntimesteps,Targets.fn(i));
    p=linspace(Xlim(1),Xlim(2),8);
    a=p(2);
    b=p(end-1);
    % choose a corner
    while(1)
        c1 = [Xlim(1)+(Xlim(2)-Xlim(1))*rand,Xlim(1)+(Xlim(2)-Xlim(1))*rand];
        c2 = [Xlim(1)+(Xlim(2)-Xlim(1))*rand,Xlim(1)+(Xlim(2)-Xlim(1))*rand];
        if norm(c1-c2)>mindist
            break
        end
    end
    
    x0=c1;
    L=c2-c1;
    magv=norm(L)/time.Tf;
    v=magv*L/norm(L);
    Targets.truth{i}(1,:)=[x0(:)',v(:)'];
    for k=2:time.Ntimesteps
        Targets.truth{i}(k,:) = Targets.f{i}(time.dt,Targets.truth{i}(k-1,:));
    end
    
    P0 = magv*diag([5^2,5^2,(0.1)^2,(0.1)^2]);
    Targets.Q{i}=diag([1^2,1^2,(0.01)^2,(0.01)^2]);
    Targets.xf{i}=zeros(time.Ntimesteps,Targets.fn(i));
    Targets.Pf{i}=zeros(time.Ntimesteps,Targets.fn(i)^2);
    Targets.Pf{i}(1,:) =  reshape(P0,1,Targets.fn(i)^2);
    Targets.xf{i}(1,:) = mvnrnd(Targets.truth{i}(1,:),P0);
end

% CT motion targets
% maxvel=1;
% minvel=0.1;
% maxTR=0.01;
% for i=9:10
%     Targets.truth{i}=zeros(time.Ntimesteps,Targets.fn(i));
%     p=linspace(Xlim(1),Xlim(2),4);
%     a=p(2);
%     b=p(end-1);
%     flg=0;
%     while(flg==0)
%         % choose a corner
%         c1=sign(randn);
%         c2=sign(randn);
%         x0=[0,0];
%         if c1<0
%             x0(1)=Xlim(1)+(a)*rand;
%         else
%             x0(1)=Xlim(2)-(a)*rand;
%         end
%         if c2<0
%             x0(2)=Xlim(1)+(a)*rand;
%         else
%             x0(2)=Xlim(2)-(a)*rand;
%         end
%         Om = -maxTR+2*maxTR*rand;
%         opc1=-1*c1;
%         opc2=-1*c2;
%         v=[opc1-c1;opc2-c2];
%         v=v/norm(v);
%         v=(minvel+(maxvel-minvel)*rand)*v;
%         Targets.truth{i}(1,:)=[x0(:)',v(:)',Om];
%         for k=2:time.Ntimesteps
%             Targets.truth{i}(k,:) = Targets.f{i}(time.dt,Targets.truth{i}(k-1,:));
%         end
%         ind1=Targets.truth{i}(:,1)>Xlim(2) | Targets.truth{i}(:,1)<Xlim(1) | Targets.truth{i}(:,2)>Xlim(2) | Targets.truth{i}(:,2)<Xlim(1);
%         if any(ind1)
%             flg=0;
%         else
%             flg=1;
%         end
%     end
%
%     P0 = diag([5^2,5^2,(0.1)^2,(0.1)^2,0.5^2]);
%     Targets.Q{i}=diag([1^2,1^2,(0.1)^2,(0.1)^2,(0.1)^2]);
%     Targets.xf{i}=zeros(time.Ntimesteps,Targets.fn(i));
%     Targets.Pf{i}=zeros(time.Ntimesteps,Targets.fn(i)^2);
%     Targets.Pf{i}(1,:) =  reshape(P0,1,Targets.fn(i)^2);
%     Targets.xf{i}(1,:) = mvnrnd(Targets.truth{i}(1,:),P0);
% end

% Now sensors - all dynamic grid based sensors
for j=1:Sensors.Nsens
    x0 = [Xlim(1)+(Xlim(2)-Xlim(1))*rand,Xlim(1)+(Xlim(2)-Xlim(1))*rand];
    xgridind = getGridindex_point(Grid,x0);
    Sensors.xfgridpos{j}=-1*ones(time.Ntimesteps,2);
    Sensors.xf{j}=-1*ones(time.Ntimesteps,2);
    
    Sensors.xf{j}(1,:) = Grid.XY(xgridind,:);
    Sensors.xfgridpos{j}(1,:)=Grid.XYgridpos(xgridind,:);
%     Sensors.xf{j}=repmat(Sensors.xf{j}(1,:),time.Ntimesteps,1);
%     Sensors.xfgridpos{j}=repmat(Sensors.xfgridpos{j}(1,:),time.Ntimesteps,1);
    
    Sensors.FOV{j}=[pi,75,pi/2]; % meters of visibility [alpha, radius, dirn]
    Sensors.h{j} = @hxy;
    Sensors.hn(j) = 2;
    Sensors.R{j} = diag([1,1]);
    Sensors.Rout{j} = diag([10,10]);
    
    
    %     [y,G,tp]= hxy(xtarg,xsenspos,FOV)
    
    Sensors.reach{j}=2; % max nodes it can jump to vert or horz from current node.
    % if reach = 0, then the sensor is static
end


figure
plot_dynamicsens(1:time.Ntimesteps,Xlim,Targets,time,Sensors,Grid,struct('Grid',1,'Sensors',1,'TargetsTruth',1,'TargetsState',0,'TargetsStateCov',0))

% keyboard

%% Running Simmulation
close all
clc
NextSensTaskTimeStep=2;
filtermethod='ut';

figure(1)
plot_dynamicsens(1:1,Xlim,Targets,time,Sensors,Grid,struct('Grid',1,'Sensors',1,'TargetsTruth',1,'TargetsState',0,'TargetsStateCov',0))

for k=2:1:time.Ntimesteps
    disp(strcat('at time step : ',num2str(k), ' : of : ',num2str(time.Ntimesteps)))
    disp('---------------------------------------------------------------------')
    tic
    % propagate
    Targets=Propagate_targets_from_Tk(Targets,time,k-1,filtermethod);
    
    %%%Generate the tasking
    if k==NextSensTaskTimeStep
        Tk=NextSensTaskTimeStep;
        TkF=NextSensTaskTimeStep+time.SensTaskHorizon;
        
        
        Sensors=DynamicSensorTask_GreedyTime_exhaust_jointcov(Targets,Sensors,Grid,time,Tk,TkF-1,filtermethod);
%         Sensors=DynamicSensorTask_GreedySensor_exhaust_jointcov_modf(Targets,Sensors,Grid,time,Tk,TkF-1,filtermethod);
%         Sensors=DynamicSensorTask_GreedySensor_exhaust_jointcov(Targets,Sensors,Grid,time,Tk,TkF-1,filtermethod);
        NextSensTaskTimeStep=min([TkF,time.Ntimesteps]);
        
    end
%     keyboard
    
    %%% Measurement update
    Targets=MeasUpdate_targets(Targets,Sensors,k,filtermethod);
    
    toc
    
    if rem(k,2)==0
        figure(1)
        clf
        plot_dynamicsens(1:k,Xlim,Targets,time,Sensors,Grid,struct('Grid',1,'Sensors',1,'TargetsTruth',1,'TargetsState',1,'TargetsStateCov',1))
        title(['time step = ',num2str(k)])
    else
        figure(2)
        clf
        plot_dynamicsens(1:k,Xlim,Targets,time,Sensors,Grid,struct('Grid',1,'Sensors',1,'TargetsTruth',1,'TargetsState',1,'TargetsStateCov',1))
        title(['time step = ',num2str(k)])
    end
%     keyboard
    pause(0.5)
    
    
    
end

save(filename)


%
% [CovMaxTrace,RMSEpos,CovFrob]=GetSatMetric(Satellites,Constants,ytruth,Constants.Ntimesteps)