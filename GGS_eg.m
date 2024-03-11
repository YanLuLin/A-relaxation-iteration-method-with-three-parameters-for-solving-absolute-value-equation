function [IT,CPU,RES]=GGS_eg(n)
eta=1.0e-6;
maxit=500;
[A,b,x_star] = eg_1(n);
[IT,CPU,RES,ERR,xk]=GGS(A,b,x_star,eta,maxit);