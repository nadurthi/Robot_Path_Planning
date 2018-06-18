function Sensors=DynamicSensorTask_GreedySensor_exhaust_jointcov(Targets,Sensors,Grid,time,Tk,TkF,method)
%%
% i need tasking from Tk to TkF both inclusive

normsqrdvec = @(X)sum(X.^2,2);


%%
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
            
%             figure(33)
%             clf
%             plot_dynamicsens(Tk:TkF,Grid.Xlim,dummyTargets,time,Sensors,Grid,struct('Grid',1,'Sensors',1,'TargetsTruth',1,'TargetsState',1,'TargetsStateCov',1))
%             hold on
%             plot( FeasibleGridPoints{k,Ntarg,nr}(:,1), FeasibleGridPoints{k,Ntarg,nr}(:,2),'m<')
%             keyboard

        end
    end
end

% keyboard

%% Find best trajectory for each sensor greedily

for nr= 1:Sensors.Nsens
    reach = Sensors.reach{nr};
    [M,counter]=gen_UAV_grid_paths(reach,Grid.dx,Ntseps);
    flg=0;
    currabspos = Sensors.xf{nr}(Tk-1,:);
    relpath=zeros(length(counter),2);
    bestpath = zeros(length(counter),2);
    bestMI=-10000000;
    cnt=0;
    while true
        cnt=cnt+1;
        if rem(cnt,1000)==0
            cnt
        end
        % get the full path of the sensor from the counter
%         counter
        for i=1:length(counter)
            relpath(i,:)=M(counter(i),2:3);
        end
        abspath=relpath2abspath_grid(currabspos,relpath,Grid);
        pruneit=prune_path(nr,abspath,FeasibleGridPoints,Grid); % pruneflg =1 means prune the shit out of it
        
        if pruneit==0 % acceptable path
            % load path into sensor path
            for k=1:Ntseps
                Sensors.xf{nr}(Tsteps(k),:)=abspath(k+1,:);
            end
            
            MI = Compute_DynamicSens_MI(Sensors,Targets,time,Tsteps,Tk);
%             keyboard
            
            if MI > bestMI
                bestMI = MI;
                bestpath = abspath;
            end
        end % prunflg
        
        if flg==1
            break
        end
        [counter,flg] = Increment_index(counter,size(M,1));
    end % exhaustive search while loop
    
    if bestMI==-10000000
        % now we have a problem, so we will move towards the target with
        % largest cov
        disp('going for path toeards large cov')
        keyboard
        bestpath = getpath2largecov(Ntseps,Tsteps,Tk,Grid,dummyTargets,Sensors,nr);
    end
    
    % save the best path
    for k=1:Ntseps
        Sensors.xf{nr}(Tsteps(k),:)=bestpath(k+1,:);
        xgridind = getGridindex_point(Grid,bestpath(k+1,:));
        Sensors.xfgridpos{nr}(Tsteps(k),:)=Grid.XYgridpos(xgridind,:);
    end
    
%     figure(33)
%     clf
%     plot_dynamicsens(Tk:TkF,Grid.Xlim,dummyTargets,time,Sensors,Grid,struct('Grid',1,'Sensors',1,'TargetsTruth',1,'TargetsState',1,'TargetsStateCov',1))
%     hold on
%     
%     keyboard
end



end