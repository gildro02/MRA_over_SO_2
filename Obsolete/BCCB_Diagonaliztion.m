N=2;
M=3;
C_0=circul(1:N);
C_1=circul(N+1:2*N);
C_2=circul(2*N+1:3*N);

C=[C_0 C_2 C_1; C_1 C_0 C_2; C_2 C_1 C_0];

W_N=(1/sqrt(N))*dftmtx(N);?
W_M=(1/sqrt(M))*dftmtx(M);

W_prod=kron(W_M,W_N);

c=reshape(C(:,1),[N,M]);
fft_c=fft2(c);
C_diagonolized=W_prod*C*W_prod';
fft_c_diag=diag(fft_c(:));
norm(C_diagonolized-fft_c_diag,"fro")