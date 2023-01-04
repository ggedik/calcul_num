function [x]=lsolve(U,b)
    [n,n] = size(U);
    x = zeros(n);
    x(1) = b(1)/U(1,1);
    for (i=2:n)
        x(i) = b(i) -U(1,1:i-1) * x(1:i-1)/U(i,i);
    end
    
endfunction
