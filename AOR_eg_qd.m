function [omega,gamma,IT,cputime,res]=AOR_eg_qd(n,maxit,omega,gamma)
tol=1.0e-6;
[A,b,x_star] = eg_1(n);
[IT,cputime,res]=AOR(A,b,x_star,omega,gamma,tol,maxit); 
