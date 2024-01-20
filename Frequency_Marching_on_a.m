a=[2,7,5,6,4i]; % k=0 to B
a_symm=[conj(flip(a(2:end))) a];
B=length(a)-1;
rho_hat=[[1/(2*pi)] 1 3 4 8]; % k=2 to B
rho_hat_extra=zeros(1,B); % B elements
rho_hat_full=[rho_hat rho_hat_extra];
rho_hat_symm=[conj(flip(rho_hat(2:end))) rho_hat];
M_1=2*pi*a_symm.*rho_hat_symm;
M_2=2*pi*diag(a_symm)*toeplitz(conj(rho_hat_full))*diag(conj(a_symm));
%sigma=0.1;
%M_2=M_2+sigma*randn(size(M_2));
k=-B:B;
k_half=0:B;
M_2_wave=M_2.*toeplitz([1./conj(M_1(k>=0)) zeros(1,B)]);

norm(M_2_wave(k>=0,k>=0)-diag(a)*toeplitz(1./conj(a))*diag(conj(a)),"fro") %should be 0

a_approx=zeros(1,B+1);
a_approx(k_half==0)=M_1(k==0);
a_approx(k_half==1)=sqrt(M_2_wave(k==1,k==1)*M_1(k==0));
a_approx_abs=sqrt(a_approx(k_half==0).*diag(M_2_wave(k>=0,k>=0)));

for n=2:B
    a_approx_raw=M_2_wave(k==n,(1<=k)&(k<n)).*...
                 flip(a_approx((1<=k_half)&(k_half<n)))./...
                 conj(a_approx((1<=k_half)&(k_half<n)));
    mean_a_approx_raw=mean(a_approx_raw);
    a_approx_phase=mean_a_approx_raw./abs(mean_a_approx_raw);
    a_approx(k_half==n)=a_approx_phase*a_approx_abs(k_half==n);
end
