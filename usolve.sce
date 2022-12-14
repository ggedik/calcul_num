function [x]=usolve(L,b)
    [n,n] = size(L);
    x = zeros(n);
    //b = L * x_ex;
    x(n)=b(n)/L(n,n);
    for(i=n-1:-1:1)
        x(i)=b(i)-L(i,(i+1):n)*x((i+1):n)/L(i,i);
    end
   // err = x_ex - x; 
endfunction
