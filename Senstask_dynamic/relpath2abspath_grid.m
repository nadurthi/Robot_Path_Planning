function abspath=relpath2abspath_grid(currabspos,relpath,Grid)
% abspath is in meters
% relpath is interms of dx from current position
n=size(relpath,1);
abspath=zeros(n+1,2);
abspath(1,:)=currabspos;
for i=1:n
    abspath(i+1,:) = abspath(i,:) + relpath(i,:);
end
