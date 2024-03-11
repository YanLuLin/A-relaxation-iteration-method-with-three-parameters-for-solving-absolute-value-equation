function [IT,CPU,RES,ERR,xk,omega,gamma,tau]=TPRI(A,b,x_star,omega,gamma,tau,eta,maxit)
xk = sparse(size(A,1),1);
D=diag(diag(A));
L=-tril(A,-1);
U=-triu(A,1);
D_L_U_inv=inv((1/omega)*D-gamma*L-tau*U);
IT = 0;
tic;
while IT<maxit
    xk =D_L_U_inv*((((1/omega)-1)*D+(1-gamma)*L+(1-tau)*U)*xk+abs(xk)+b);
    RES=norm(b+abs(xk)-A*xk)/norm(b);
    if RES<eta          
        IT = IT+1;
        break
    end
    IT = IT+1;
end
CPU=toc;
ERR=norm(xk-x_star);