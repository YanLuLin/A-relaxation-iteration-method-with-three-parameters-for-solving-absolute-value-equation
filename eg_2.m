function [A,b,x_star] = eg_2(m,mu)
n=m^2;
S=sparse(4*speye(m)-diag(ones(m-1,1),-1)-diag(ones(m-1,1),1));
A_hat=kron(-diag(ones(1,m-1),-1)-diag(ones(1,m-1),1),speye(m))+kron(diag(ones(1,m)),S);  
A=sparse(A_hat+mu*speye(n));
x_star=zeros(n,1)+(-1).^(1:n)'; 
b=A*x_star-abs(x_star);