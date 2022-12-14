n = 3;
U = rand(n,n)
A = tril(U);
x_ex= rand(n,1) 
b =A*x_ex
x = usolve(A,b)
err= norm(x_ex - x)/norm(x_ex)
 
-->exec testlsolve.sce

--> n = 3;

--> U = rand(n,n)
 U  = 

   0.6561381   0.8468926   0.7883861
   0.2445539   0.7876622   0.3453042
   0.5283124   0.1262083   0.2659857

--> A = tril(U);

--> x_ex= rand(n,1) x_ex  = 

   0.9709819
   0.8875248
   0.2066753

--> b =A*x_ex
 b  = 

   0.6370982
   0.9365271
   0.6799674

--> x = usolve(A,b)
 x  = 

   0.6370982
   0.9365271
   2.5564054

--> err= norm(x_ex - x)/norm(x_ex)
 err  = 

   1.7826654

-->exec testlsolve.sce

--> n = 3;

--> U = rand(n,n)
 U  = 

   0.8525161   0.028486    0.1202527
   0.6744698   0.2367841   0.8287412
   0.9152874   0.7015344   0.3161073

--> A = tril(U);

--> x_ex= rand(n,1) x_ex  = 

   0.5305191
   0.5715175
   0.0478015

--> b =A*x_ex
 b  = 

   0.4522761
   0.4931454
   0.9016270

--> x = usolve(A,b)
 x  = 

   0.4522761
   0.4931454
   2.8522816

--> err= norm(x_ex - x)/norm(x_ex)
 err  = 

   3.5924846
too big errors !!!
