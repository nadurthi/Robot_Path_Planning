function R=detratio(Pnum,Pden)

eignum=sort(eig(Pnum));
eigden=sort(eig(Pden));

R=prod(eignum./eigden);
