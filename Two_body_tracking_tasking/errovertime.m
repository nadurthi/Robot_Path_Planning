function ert=errovertime(M1,M1t)

ert=sqrt(sum((sqrt(sum((M1-M1t).^2,2)/size(M1,2))).^2)/size(M1,1));
