function [x] = usolve(U,b)
    x(1) = b(1)/U(1,1);
for i=2:n
    x(i)=(b(i)-U(i,1:(i-1))*x(1:(i-1)))/U(i,i);
end
endfunction
