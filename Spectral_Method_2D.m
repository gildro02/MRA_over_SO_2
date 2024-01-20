B=50;
Q=5;
M=2*B+1;
%a=exp(1:s).'.*(1:s).';
%a=a./abs(a);
rho_hat=[1/(2*pi); (1:2*B).'];
C=BCCB(repmat(rho_hat.',Q,1));
T_blocks=cell(M,M);
T_sub_matrices=cell(M,1);
for n=1:M
    T_sub_matrices{n}=exp(-((n-1)^(1/8))).*ones(Q,Q).*((n-1)<=M/1.5);
end
for n=1:M
    for m=1:M
        T_blocks{n,m}=T_sub_matrices{abs(n-m)+1};
    end
end
T=cell2mat(T_blocks);


W_M=(1/sqrt(M))*dftmtx(M);
W_Q=(1/sqrt(Q))*dftmtx(Q);
W=kron(W_M,W_Q);

a=exp(7i*(1:Q*M)).';
%a=(1:Q*M).';
a_tilde=a./abs(a);
A=BCCB(ifft2(mat(a_tilde,Q,B)));
P_a=abs(a).^2;
v_p=1./sqrt(P_a);
M1_true=a.*C(:,1); %only for use of C, not T.
M2_true=diag(a)*C*diag(a)';
M2_approx=diag(a)*T*diag(a)';
M2_mathcal=W'*diag(v_p)*M2_true*diag(v_p)*W;

[V,D]=eig(M2_mathcal,"vector");
[~,u_index]=max(min(abs(D-D.')+diag(inf(1,Q*M)),[],2)); %unique eigenvector (index)

v_kappa=V(:,u_index); %unique eigenvector
a_est_tilde=vec(fft2(mat(v_kappa,Q,B)));
a_est = sqrt(P_a).* a_est_tilde;
a_est=a_est.*a_tilde(1)/a_est_tilde(1);
[err_squared,l]=Circ_Error_Continuous_2D(a_est,a,Q,B)