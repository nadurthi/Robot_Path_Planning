

A=zeros(1,5);
cnt=1;
tic
while(1)
    [A,flg] = Increment_index(A,15);
    if flg==1
        break
    end
   cnt=cnt+1;
   if rem(cnt,10000)==0
       cnt
   end
end

toc