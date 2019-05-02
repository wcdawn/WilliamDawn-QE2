clear
close all
rng(1);

N = 3000;

A = gallery('condex',N);
b = rand(N,1);

fprintf('cond_two = %.6e\n',cond(A,2));
fprintf('cond_inf = %.6e\n',cond(A,inf));
tic;
x_matlab = A\b;
toc;

tol = 1d-13;
maxit = 100;
matvec = @(x)(A*x);
tic;
[x] = bicgstab(b, matvec, tol, maxit);
toc;

fprintf('res_matlab = %.6e\n',norm(b-A*x_matlab,2))
fprintf('res        = %.6e\n',norm(b-A*x,2))

err = x-x_matlab;
fprintf('err_norm_two = %.6e\n',norm(err,2));
fprintf('err_norm_inf = %.6e\n',norm(err,inf));

figure
plot(x_matlab);
hold on
plot(x);
legend({'x_{matlab}','x'});
