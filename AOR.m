 function [IT,CPU,RES,ERR,xk]=AOR(A,b,x_star,omega,gamma,eta,maxit)
xk = sparse(size(A,1),1);
D=diag(diag(A));
L=-tril(A,-1);
U=-triu(A,1);
D_rL_inv=inv(D-gamma*L);
IT = 0;
tic;
while IT<maxit
    xk =D_rL_inv*(((1-omega)*D+(omega-gamma)*L+omega*U)*xk+omega*abs(xk)+omega*b);
    RES=norm(b+abs(xk)-A*xk)/norm(b);
    if RES<eta          
        IT = IT+1;
        break
    end
    IT = IT+1;
end
CPU=toc;
ERR=norm(xk-x_star);