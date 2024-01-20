Kappa=100;
B=10;
seed_rho_coeff=(B/2):-1:1;
[rho_coeff,rho]=Generate_Rho_Squared(seed_rho_coeff);
rho_coeff_compressed=rho_coeff(ceil(end/2):end);
coeff_compressed=Kappa:-1:1;
lambda=1e-9;
pad=B;
coeff_padded=[coeff_compressed lambda*ones(1,pad)];

T_bad=toeplitz([conj(rho_coeff_compressed) zeros(1,Kappa+pad-B-1)]);
T_good=toeplitz([conj(rho_coeff_compressed) zeros(1,Kappa+pad-2*B-1) flip(rho_coeff_compressed(2:end))]);

%c=[rho_coeff_compressed zeros(1,Kappa+pad-B-1)];
%T=toeplitz(c,[c(1) flip(conj(c(2:end)))]);


M_2_good=diag(coeff_padded)*T_good*diag(conj(coeff_padded));
M_2_bad=diag(coeff_padded)*T_bad*diag(conj(coeff_padded));

%M_2_test=(coeff_compressed.'*conj(coeff_compressed)).*T;

W=dftmtx(Kappa+pad);
M_2_classic_good=W*M_2_good*W';
M_2_classic_bad=W*M_2_bad*W';

v_p=1./abs(coeff_padded);
Q=W*diag(v_p)*W';

M_2_classic_wave_good=Q*M_2_classic_good*Q';
M_2_classic_wave_bad=Q*M_2_classic_bad*Q';

[V_good,D_good]=eig(M_2_classic_wave_good,"vector");
[V_bad,D_bad]=eig(M_2_classic_wave_bad,"vector");

%pick the eigenvector corresponding to the "most isolated" eigenvalue:
[~,u_index_good]=max(min(abs(D_good-D_good.')+diag(inf(1,Kappa+pad)),[],2)); 
[~,u_index_bad]=max(min(abs(D_bad-D_bad.')+diag(inf(1,Kappa+pad)),[],2)); 

coeff_recovered_good=fft(V_good(:,u_index_good))./v_p.'; %abs to get rid of global phase, no reason to assume realness
coeff_recovered_bad=fft(V_bad(:,u_index_bad))./v_p.';

abs_coeff_recovered_good=abs(coeff_recovered_good);
abs_coeff_recovered_bad=abs(coeff_recovered_bad);

%coeff_recovered_true_good=coeff_recovered./(coeff_recovered(1)/coeff_padded(1));
%coeff_recovered_true_bad=coeff_recovered./(coeff_recovered(1)/coeff_padded(1));

coeff_recovered_true_good=abs(coeff_recovered_good);
coeff_recovered_true_bad=abs(coeff_recovered_bad);

error_good=norm(coeff_recovered_true_good(1:Kappa)-coeff_padded(1:Kappa).')./norm(coeff_padded)
error_bad=norm(coeff_recovered_true_bad(1:Kappa)-coeff_padded(1:Kappa).')./norm(coeff_padded)
