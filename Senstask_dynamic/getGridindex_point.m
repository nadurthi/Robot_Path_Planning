function ind = getGridindex_point(Grid,x)
[m,ind]=min(sum((Grid.XY-repmat(x(:)',size(Grid.XY,1),1)).^2,2));
