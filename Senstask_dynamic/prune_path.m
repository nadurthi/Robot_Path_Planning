function pruneit=prune_path(ns,abspath,FeasibleGridPoints,Grid)
% ns is the sensor for which we want to check pruning
% remember abspath starts from Tk-1
thresh =0.01;
threshsqred = thresh^2;
%% first check if path is within grid
if any(abspath(:,1)>Grid.Xlim(2)) | any(abspath(:,1)<Grid.Xlim(1)) | any(abspath(:,2)>Grid.Xlim(2)) | any(abspath(:,2)<Grid.Xlim(1))
    pruneit=1;
    return
end


%% then check if paths + circular FOV intersectt the covs

% FeasibleGridPoints=cell(Ntseps,Targets.Ntargs,Sensors.Nsens)
[Ntsteps,Ntargs,Nsens] = size(FeasibleGridPoints);
pruneit=1;
for k=1:Ntsteps
    for Ntarg=1:Ntargs
        N = size(FeasibleGridPoints{k,Ntarg,ns},1);
        
        
        if any(ismember(abspath(k+1,:),FeasibleGridPoints{k,Ntarg,ns})==1)
           pruneit=0;
           break;
        end
    end
    if pruneit==0
        break
    end
end



















