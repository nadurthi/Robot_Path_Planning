function abspath = getpath2largecov(Ntsteps,Tsteps,Tk,Grid,dummyTargets,Sensors,nr)
% Just generate a path towards the target with largest cov 
% Tsteps = Tk:TkF;
% greedily pick steps that get you closer to the target with largest cov


abspath=zeros(Ntsteps+1,2);
abspath(1,:)=Sensors.xf{nr}(Tk-1,:);

maxtargtime=[];
maxcov=-10000;
for k=1:Ntsteps
    for Ntarg=1:dummyTargets.Ntargs
        %                 mk=dummyTargets.xf{Ntarg}(Tsteps(k),:)';
        Pk=reshape(dummyTargets.Pf{Ntarg}(Tsteps(k),:),dummyTargets.fn(Ntarg),dummyTargets.fn(Ntarg));
        if max(trace(Pk))>maxcov
            maxtargtime = [k,Ntarg];
            maxcov=max(trace(Pk));
        end
    end
end
mk=dummyTargets.xf{maxtargtime(2)}(Tsteps(maxtargtime(1)),:)';
% bestline = mk(1:2)'- Sensors.xf{nr}(Tk-1,:);

reach = Sensors.reach{nr}; % number of dx units it can move to left and right: it is square
[M,~]=gen_UAV_grid_paths(reach,Grid.dx,Ntsteps);
N=size(M,1);
for k=1:Ntsteps
    currpos = abspath(k,:);
    bestpos=[];
    dis=1000000000000;
    for i=1:N
        mpos = currpos+M(i,2:3);
        dist2targ = norm(mk(1:2)'-mpos);
        if dist2targ<dis
            dis= dist2targ;
            bestpos = mpos;
        end
    end
    abspath(k+1,:)=bestpos;
end