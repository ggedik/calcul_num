x = [1,2,3,4];
y = [1;2;3;4];
z = x*y;
size(x);
size(y);
y=norm(x,2);
A= [[1,2,3];[3,4,5];[6,7,8];[9,10,11]];
A= [[1,2,3];[3,4,5];[6,7,8]];
cond(A);
A = rand(3,3); //matrice randomly rempli 3x3
x = rand(1,3);
x = rand(3,1);
x_ex = rand(3,1);
b = A*x_ex;
x_app = A\b;
err_av = norm(x_ex-x_app)/norm(x_ex);
err_ar = norm(b-A*x_app)/norm(A)*norm(x_app);
 
A = rand(1000,1000); //matrice randomly rempli 1k
x_ex = rand(1000,1);
b = A*x_ex;
x_app = A\b;
err_av = norm(x_ex-x_app)/norm(x_ex);
err_ar = norm(b-A*x_app)/norm(A)*norm(x_app);

A = rand(100,100); //matrice randomly rempli 10k, takes too much time 
x_ex = rand(100,1);
b = A*x_ex;
x_app = A\b;
err_av = norm(x_ex-x_app)/norm(x_ex);
err_ar = norm(b-A*x_app)/norm(A)*norm(x_app);

A=rand(5,5);
tril(A);
