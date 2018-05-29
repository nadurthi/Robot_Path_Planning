function [A,flg] = Increment_index(A,max_ind)
% right most or last index in A is first incremented

if min(A)==max_ind
    flg=1;
    return
end


n=length(A);
flg=0;
while(1)
    a = A(n)+1;
    if a<=max_ind
        A(n)=a;
        break;
    else
        A(n)=0;
        n=n-1;
        
    end
    
    if n<=0
        flg=1;
        break
    end
    
    
    
end