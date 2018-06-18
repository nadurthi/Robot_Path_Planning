function Sensors=DynamicSensorTask_GreedyTime_exhaust_jointcov(Targets,Sensors,Grid,time,Tk,TkF,method)
%%
% i need tasking from Tk to TkF both inclusive
% greedy time



normsqrdvec = @(X)sum(X.^2,2);


%%

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


Tsteps = Tk:TkF;
Ntseps = length(Tsteps);
for nr= 1:Sensors.Nsens
    Sensors.xf{nr}(Tsteps,:)=-1;
end

%% First propagate the targets to get the covariance ellipsoids
% used to prune the trajectories
dummyTargets = Targets;
for k=Tk:TkF-1
    dummyTargets=Propagate_targets_from_Tk(dummyTargets,time,k,method);
end
FeasibleGridPoints=cell(Ntseps,Targets.Ntargs,Sensors.Nsens);

Nmcpts=1000;
Ngrid = size(Grid.XY,1);

     
% Now get the feasible points of visibility
for k=1:Ntseps
    
    for Ntarg=1:dummyTargets.Ntargs
        mk=dummyTargets.xf{Ntarg}(Tsteps(k),:)';
        Pk=reshape(dummyTargets.Pf{Ntarg}(Tsteps(k),:),dummyTargets.fn(Ntarg),dummyTargets.fn(Ntarg));
        
        mk=mk(1:2);
        Pk=Pk(1:2,1:2);
        X=mvnrnd(mk',Pk,Nmcpts);
        eigmax = max(eig(Pk));
        for nr= 1:Sensors.Nsens
            
            FeasibleGridPoints{k,Ntarg,nr}=[];
            Rmax = Sensors.FOV{nr}(2);
            r = eigmax+Rmax;
            potential_pts = Grid.XY( normsqrdvec(Grid.XY-repmat(mk',Ngrid,1))<=r^2,:);
            for i=1:size(potential_pts,1)
                S=normsqrdvec(X-repmat(potential_pts(i,:),Nmcpts,1))<=Rmax^2;
                if sum(S)/Nmcpts >= 0.25
                    FeasibleGridPoints{k,Ntarg,nr} = [FeasibleGridPoints{k,Ntarg,nr};potential_pts(i,:)] ;
                end
            end
            


        end
    end
end

% keyboard
%% get the reachability sets for each sensor and at each time step.
ReachabilitySetPoints=cell(Sensors.Nsens,1+Ntseps);

for nr= 1:Sensors.Nsens
    reach = Sensors.reach{nr};
    currabspos = Sensors.xf{nr}(Tk-1,:);
    ReachabilitySetPoints{nr,1}=urrabspos;
    [M,counter]=gen_UAV_grid_paths(reach,Grid.dx,Ntseps);
    for k=1:Ntseps
        N=size(ReachabilitySetPoints{nr,k},1);
        ReachabilitySetPoints{nr,k+1}=[];
        for j=1:N 
            currabspos = ReachabilitySetPoints{nr,k}(j,:);
            for i=1:size(M,1)
                ReachabilitySetPoints{nr,k+1} = [ReachabilitySetPoints{nr,k+1};currabspos+M(i,2:3)];
            end
        end
        ReachabilitySetPoints{nr,k+1}=unique(ReachabilitySetPoints{nr,k+1},'rows');
        
        figure(33)
        clf
        plot_dynamicsens(Tk:TkF,Grid.Xlim,dummyTargets,time,Sensors,Grid,struct('Grid',1,'Sensors',1,'TargetsTruth',1,'TargetsState',1,'TargetsStateCov',1))
        hold on
        plot( ReachabilitySetPoints{1,k+1}(:,1), ReachabilitySetPoints{nr,k+1}(:,2),'m<')
        plot( ReachabilitySetPoints{2,k+1}(:,1), ReachabilitySetPoints{nr,k+1}(:,2),'cs','MarkerSize',6)
        keyboard

    end
end

%% Greedy in time
% first greedily add 
for k=Ntseps:1
    max_ind_vec=zeros(1,Sensors.Nsens);
    S=cell(1,Sensors.Nsens);
    for nr= 1:Sensors.Nsens
        if k==Ntseps
            S{nr}=ReachabilitySetPoints{nr,k+1};
        else
            reach = Sensors.reach{nr};
            [M,~]=gen_UAV_grid_paths(reach,Grid.dx,Ntseps);
            currpos=Sensors.xf{nr}(Tsteps(k+1),:);
            for m=1:size(M,1)
               M(i,:)=currpos+M(i,:); 
            end
            ind= M(:,1)>Grid.Xlim(2) | M(:,1)<Grid.Xlim(1) | M(:,2)>Grid.Xlim(2) | M(:,2)<Grid.Xlim(1);
            M(ind,:)=[];
            S{nr}=intersect(ReachabilitySetPoints{nr,k+1},M,'rows');
        end
        max_ind_vec(nr)=size(S{nr},1);
    end
    A=ones(1,Sensors.Nsens);
    flg=1;
    bestMI=0;
    bestPos=cell(1,Sensors.Nsens);
    
    figure(33)
    clf
    plot_dynamicsens(Tk:TkF,Grid.Xlim,dummyTargets,time,Sensors,Grid,struct('Grid',1,'Sensors',1,'TargetsTruth',1,'TargetsState',1,'TargetsStateCov',1))
    hold on
    plot( S{1}(:,1), S{nr}(:,2),'m<')
    plot( S{2}(:,1), S{nr}(:,2),'cs','MarkerSize',6)
    keyboard
        
    while true
        for nr= 1:Sensors.Nsens
            Sensors.xf{nr}(Tsteps(k),:) = S{nr}(A(nr),:);
        end
        MI = Compute_DynamicSens_MI(Sensors,Targets,time,Tsteps,Tk);
        if MI>bestMI
            bestMI=MI;
            for nr= 1:Sensors.Nsens
                bestPos{nr}=Sensors.xf{nr}(Tsteps(k),:);
            end
        end
        if flg==1
            break
        end
        [A,flg] = Increment_index_indivmax(A,max_ind_vec);
    end
    for nr= 1:Sensors.Nsens
        Sensors.xf{nr}(Tsteps(k),:)=bestPos{nr};
    end
end



end