 function [IT,CPU,RES,ERR,xk]=GAOR(A,b,x_star,u,omega,gamma,eta,maxit)
xk = sparse(size(A,1),1);
O=u*eye(size(A));
D=diag(diag(A));
L=-tril(A,-1);
U=-triu(A,1);
A_inv=inv(A);
D_rL_inv=inv(O+D-gamma*L);
IT = 0;
tic;
while IT<maxit
    yk=A_inv*(abs(xk)+b);
    xk =D_rL_inv*((O+(1-omega)*D+(omega-gamma)*L+omega*U)*xk+omega*abs(yk)+omega*b);
    RES=norm(b+abs(xk)-A*xk)/norm(b);
    if RES<eta         
        IT = IT+1;
        break
    end
    IT = IT+1;
end
CPU=toc;
ERR=norm(xk-x_star);