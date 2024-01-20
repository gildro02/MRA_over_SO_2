Kappa=100;
B=10;
seed_rho_coeff=(B/2):-1:1;
[rho_coeff,rho]=Generate_Rho_Squared(seed_rho_coeff);
rho_coeff_compressed=rho_coeff(ceil(end/2):end);
coeff_compressed=Kappa:-1:1;
%T=toeplitz([conj(rho_coeff_compressed) zeros(1,Kappa-B-1)]);
c=[rho_coeff_compressed zeros(1,Kappa-B-1)];
T=toeplitz(c,[c(1) flip(conj(c(2:end)))]);
M_2=diag(coeff_compressed)*T*diag(conj(coeff_compressed));
%M_2_test=(coeff_compressed.'*conj(coeff_compressed)).*T;

W=dftmtx(Kappa);
M_2_classic=W*M_2*W';
v_p=1./abs(coeff_compressed);
Q=W*diag(v_p)*W';
M_2_classic_wave=Q*M_2_classic*Q';
[V,D]=eig(M_2_classic_wave,"vector");

%pick the eigenvector corresponding to the "most isolated" eigenvalue:
[~,u_index]=max(min(abs(D-D.')+diag(inf(1,Kappa)),[],2)); 
%u_index=96;
coeff_recovered=abs(fft(V(:,u_index))./v_p.'); %abs to get rid of global phase, no reason to assume realness

error=Circular_Error(coeff_recovered.',coeff_compressed)./norm(coeff_compressed);