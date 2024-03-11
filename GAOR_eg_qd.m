function [omega,gamma,IT,cputime,res]=GAOR_eg_qd(u,n,maxit,omega,gamma)
tol=1.0e-6;
[A,b,x_star] = eg_1(n);
[IT,cputime,res]=GAOR(A,b,x_star,u,omega,gamma,tol,maxit); 
