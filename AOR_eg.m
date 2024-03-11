function [opt_omega,opt_gamma,opt_IT,opt_cputime,res]=AOR_eg(n,maxit,omega_0,omega_N,delta_omega,gamma_0,gamma_N,delta_gamma)
tol=1.0e-6;
[A,b,x_star] = eg_1(n);
N_omega=floor((omega_N-omega_0)/delta_omega);       
N_gamma=floor((gamma_N-gamma_0)/delta_gamma);        
for j=1:N_omega
    omega(j)=omega_0+j*delta_omega;
    for k=1:N_gamma
        gamma(j,k)=gamma_0+k*delta_gamma;
        [IT(j,k),cputime(j,k),res(j,k)]=AOR(A,b,x_star,omega(j),gamma(j,k),tol,maxit); 
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[opt_IT,location]=find_min(IT);
for q=1:size(location,1)
        L1=location(q,1);
        L2=location(q,2);
        min_IT_omega(q)=omega(L1);
        min_IT_gamma(q)=gamma(L1,L2);
        min_IT_cputime(q)=cputime(L1,L2);
        min_IT_res(q)=res(L1,L2);
end
min_IT_para=[min_IT_omega',min_IT_gamma'];
[opt_cputime,k]=min(min_IT_cputime);
opt_omega=min_IT_omega(k);
opt_gamma=min_IT_gamma(k);
res=min_IT_res(k);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%