function [opt_omega,opt_IT,opt_cputime,res]=SOR_eg(n,maxit,omega_0,omega_N,delta_omega)
tol=1.0e-6;
[A,b,x_star] = eg_1(n);
N_omega=floor((omega_N-omega_0)/delta_omega);        
for j=1:N_omega
    omega(j)=omega_0+j*delta_omega;
    [IT(j),cputime(j),res(j)]=SOR(A,b,x_star,omega(j),tol,maxit);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(omega,IT)
opt_IT=min(IT);
location=find(IT==opt_IT);
for q=1:size(location)
        L=location(q);
        min_IT_omega(q)=omega(L);
        min_IT_cputime(q)=cputime(L);
        min_IT_res(q)=res(L);
end
min_IT_omega
[opt_cputime,k]=min(min_IT_cputime);
opt_omega=min_IT_omega(k);
res=min_IT_res(k);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [opt_cputime,k]=min(cputime);
% opt_IT=IT(k);
% opt_omega=omega(k);