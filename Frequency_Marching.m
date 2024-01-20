a=[1,7,5,6,4i]; % k=0 to B
a_symm=[conj(flip(a(2:end))) a];
rho_hat=[[1] 1 3 4]; % k=2 to B
rho_hat_extra=[5 6 7]; % B elements
rho_hat_full=[rho_hat rho_hat_extra];
rho_hat_symm=[conj(flip(rho_hat(2:end))) rho_hat];
rho_hat_full_symm=[conj(flip(rho_hat(2:end))) rho_hat];
B=length(a)-1;
M_1=a_symm.*rho_hat_symm;
M_2=diag(a_symm)*toeplitz(conj(rho_hat_full))*diag(conj(a_symm));

a_approx=zeros(1,B+1);
a_approx(1)=M_1(B+1);
rho_hat_approx=zeros(1,B+1);
rho_hat_approx(1)=1;
rho_hat_1=1; %assumption, we have this freedom in the first coeff
diag_minus_1=diag(M_2,-1);
for b=1:B
    a_approx(b+1)=diag_minus_1(B+b)/(a_approx(b)'*rho_hat_approx(1));
end
a_approx_symm=[conj(flip(a_approx(2:end))) a_approx];
rho_hat_approx_full_symm=M_1./a_approx_symm;
