function [IT,CPU,RES,ERR,xk]=GGS(A,b,x_star,eta,maxit)
xk = sparse(size(A,1),1);
D=diag(diag(A));
L=-tril(A,-1);
U=-triu(A,1);
A_inv=inv(A);
D_L_inv=inv(D-L);
IT = 0;
tic;
while IT<maxit
    yk=A_inv*(abs(xk)+b);
    xk =D_L_inv*(U*xk+abs(yk)+b);
    RES=norm(b+abs(xk)-A*xk)/norm(b);
    if RES<eta          
        IT = IT+1;
        break
    end
    IT = IT+1;
end
CPU=toc;
ERR=norm(xk-x_star);