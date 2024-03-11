 function [IT,CPU,RES,ERR,xk]=SOR(A,b,x_star,omega,eta,maxit)
xk = sparse(size(A,1),1);
D=diag(diag(A));
L=-tril(A,-1);
U=-triu(A,1);
wD_L_inv=inv((1/omega)*D-L);
IT = 0;
tic;
while IT<maxit
    xk =wD_L_inv*((((1/omega)-1)*D+U)*xk+abs(xk)+b);
    RES=norm(b+abs(xk)-A*xk)/norm(b);
    if RES<eta         
        IT = IT+1;
        break
    end
    IT = IT+1;
end
CPU=toc;
ERR=norm(xk-x_star);