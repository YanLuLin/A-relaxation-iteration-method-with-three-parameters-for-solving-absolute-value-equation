function [opt_omega,opt_gamma,opt_tau,opt_IT,opt_cputime,res]=TPRI_eg(n,maxit,omega_0,omega_N,delta_omega,gamma_0,gamma_N,delta_gamma,tau_0,tau_N,delta_tau)
tol=1.0e-6;
[A,b,x_star] = eg_1(n);
N_omega=floor((omega_N-omega_0)/delta_omega);        
N_gamma=floor((gamma_N-gamma_0)/delta_gamma);        
N_tau=floor((tau_N-tau_0)/delta_tau);                
for j=1:N_omega
    omega(j)=omega_0+j*delta_omega;
    for k=1:N_gamma
        gamma(j,k)=gamma_0+k*delta_gamma;
        for i=1:N_tau
            tau(j,k,i)=tau_0+i*delta_tau;
            [IT(j,k,i),cputime(j,k,i),res(j,k,i)]=TPRI(A,b,x_star,omega(j),gamma(j,k),tau(j,k,i),tol,maxit);
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opt_cputime=min(cputime(:));
s=size(cputime);
Lin=find(cputime<=opt_cputime);
Lin=Lin(1);
[j,k,i]=ind2sub(s,Lin);
opt_omega=omega(j);
opt_gamma=gamma(j,k);
opt_tau=tau(j,k,i);
opt_IT=IT(j,k,i);
res=res(j,k,i);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
