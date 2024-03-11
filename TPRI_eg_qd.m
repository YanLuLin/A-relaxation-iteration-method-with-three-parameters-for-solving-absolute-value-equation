function [omega,gamma,tau,IT,CPU,RES]=TPRI_eg_qd(n,maxit,omega,gamma,tau)
eta=1.0e-6;
[A,b,x_star] = eg_1(n);
[IT,CPU,RES,ERR,xk,omega,gamma,tau]=TPRI(A,b,x_star,omega,gamma,tau,eta,maxit); 
