function P=rndCov(n)
P=randn(n,n);
P=(P*P').^2;
end