function [A,flg] = Increment_binary_index(A)
% right most or last index in A is first incremented

if min(A)==1
    flg=1;
    return
end


n=length(A);
flg=0;
while(1)
    a = A(n)+1;
    if a<=1
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

if min(A)==1
    flg=1;
    return
end