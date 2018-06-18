function [M,counter]=gen_UAV_grid_paths(reach,dx,Nt)
% reachability: number of cartesian consecutive nodes, it is rectangle
% for Nt time steps
% 

L= Nt*reach;

r=2*reach+1;
M=zeros(r^2,3);
i=1;
j=1;
k=1;
while 1
    M(k,1) = k;
    M(k,2) = (j-1-reach)*dx;
    M(k,3) = -(i-1-reach)*dx;
    k=k+1;
    
    j=j+1;
    if j>r
        j=1;
        i=i+1;
    end
    if i>r
        break
    end
end

counter = ones(1,Nt);
% Now call
%[counter,flg] = Increment_index(counter,size(M,1));



