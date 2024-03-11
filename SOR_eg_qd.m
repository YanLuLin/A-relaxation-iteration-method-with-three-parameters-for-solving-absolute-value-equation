function [omega,IT,CPU,RES]=SOR_eg_qd(n,maxit,omega)
eta=1.0e-6;
[A,b,x_star] = eg_1(n);
[IT,CPU,RES,ERR,xk]=SOR(A,b,x_star,omega,eta,maxit);
