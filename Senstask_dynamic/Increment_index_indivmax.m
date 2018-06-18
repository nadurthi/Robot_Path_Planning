function [A,flg] = Increment_index_indivmax(A,max_ind_vec)
% right most or last index in A is first incremented

if all(max_ind_vec-A==0)
    flg=1;
    return
end


n=length(A);
flg=0;
while(1)
    a = A(n)+1;
    if a<=max_ind_vec(n)
        A(n)=a;
        break;
    else
        A(n)=1;
        n=n-1;
        
    end
    
    if n<=0
        flg=1;
        break
    end
    
    
    
end