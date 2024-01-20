n=1000;
x=((n:-1:1).');
v_p=1./sqrt(n)*1./abs(fft(x));
k=470;
%{
rho=[k:-1:1 zeros(1,n-k)].';
rho=rho./sum(rho);
%B=10;
%rho=ifft(fft(rho).*[ones(B,1); zeros(n-B,1)]);
M1=circul(rho)*x;
M2=circul(x)*diag(rho)*circul(x).';
[x_est,rho_est]=MRA_Single_Ring(v_p,M1,M2,0);
Circular_Error(x_est,x)
Circular_Error(rho_est,rho)
%}

rho_fft_seed=sqrt((k/2:-1:1).');
%{
rho_fft_seed_real=[rho_fft_seed; zeros()]
rho_fft=cconv([rho_fft_seed; zeros(n-length(rho_fft_seed),1)],...
        [rho_fft_seed; zeros(n-length(rho_fft_seed),1)],n);
rho=ifft(rho_fft);
%}
rho_coeff=Generate_Rho_Squared(rho_fft_seed);
rho_coeff_half=rho_coeff((end+1)/2:end);
T=toeplitz([conj(rho_coeff((end+1)/2:end)); zeros(n-length(rho_coeff_half),1)]);
W=(1/sqrt(n))*dftmtx(n);
M2_T=W'*diag(fft(x))*T*diag(conj(fft(x)))*W;
M1_T=circul(x)*ifft(rho_coeff_half,n);
[x_T,rho_T]=MRA_Single_Ring(v_p,M1_T,M2_T,0);
