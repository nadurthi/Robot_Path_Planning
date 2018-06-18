function plot_dynamicsens(Tfs,Xlim,Targets,time,Sensors,Grid,pltflags)
% struct('Grid',1,'Sensors',1,'TargetsTruth',1,'TargetsState',1,'TargetsStateCov',1)
% Tfs is time steps
% Sensors.FOV{j}=[pi,50,pi/2]; % meters of visibility [alpha, radius, dirn]

hold on

% plot grid
if pltflags.Grid==1
    plot(Grid.XY(:,1),Grid.XY(:,2),'k.')
end

% plot UAVS
if pltflags.Sensors==1
    for i=1:Sensors.Nsens
        ind=Sensors.xf{i}(Tfs,1)>=0;
        SS=Tfs(ind);
        plot(Sensors.xf{i}(SS,1),Sensors.xf{i}(SS,2),'k');
        plot(Sensors.xf{i}(SS(end),1),Sensors.xf{i}(SS(end),2),'ks','MarkerSize',6);
        plot_circlar_sensor(Sensors.xf{i}(SS(end),1),Sensors.xf{i}(SS(end),2),Sensors.FOV{i}(2),'g')
    end
end


% plot targets
if pltflags.TargetsTruth==1
    for i=1:Targets.Ntargs
        plot(Targets.truth{i}(Tfs,1),Targets.truth{i}(Tfs,2),'b--');
        plot(Targets.truth{i}(Tfs(end),1),Targets.truth{i}(Tfs(end),2),'bo','MarkerSize',6);
    end
end


% plot targets state
if pltflags.TargetsState==1
    for i=1:Targets.Ntargs
        plot(Targets.xf{i}(Tfs,1),Targets.xf{i}(Tfs,2),'r--');
        plot(Targets.xf{i}(Tfs(end),1),Targets.xf{i}(Tfs(end),2),'r*','MarkerSize',6);
    end
end

% plot target cov
if pltflags.TargetsStateCov==1
    for i=1:Targets.Ntargs
        Px=reshape(Targets.Pf{i}(Tfs(end),:),Targets.fn(i),Targets.fn(i));
        plot_targ_ellipse(Targets.xf{i}(Tfs(end),1:2),Px(1:2,1:2));
        plot(Targets.xf{i}(Tfs(end),1),Targets.xf{i}(Tfs(end),2),'r*','MarkerSize',6);
    end
end
axis([Xlim(1)-100,Xlim(2)+100,Xlim(1)-100,Xlim(2)+100])
% axis equal
axis square
hold off
