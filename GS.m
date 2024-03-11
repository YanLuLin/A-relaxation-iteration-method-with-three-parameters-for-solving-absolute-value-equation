function [IT,CPU,RES,ERR,xk]=GS(A,b,x_star,eta,maxit)
xk = sparse(size(A,1),1);
D=diag(diag(A));
L=-tril(A,-1);
U=-triu(A,1);
D_L_inv=inv(D-L);
IT = 0;
tic;
while IT<maxit
    xk =D_L_inv*(U*xk+abs(xk)+b);
    RES=norm(b+abs(xk)-A*xk)/norm(b);
    if RES<eta          
        IT = IT+1;
        break
    end
    IT = IT+1;
end
CPU=toc;
ERR=norm(xk-x_star);