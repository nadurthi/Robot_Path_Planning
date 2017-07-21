function [X,iter]=getitertensorprod(possSens,Nrads,Nr,iter)
% possSens is the possibilities for the Nrth radars
% there are total of Nrads
% iter is the vector [3,2,4,4] .. of length Nrads,

Np=size(possSens,1);
if iter(Nr) > Np
    iter(Nr)=1;
end