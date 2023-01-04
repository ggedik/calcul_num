n=3;
U = rand(n,n)
A = triu(U);
x_ex= rand(n,1) 
b =A*x_ex
x = usolve(A,b)
err= norm(x_ex - x)/norm(x_ex)
