c=[4 2 1].';
Q=[3 2 1].';
s=sum(Q);
M=length(Q);
a=exp(1:s).'.*(1:s).';
a=a./abs(a);
%{
C=[a(1)*ones(Q(1),Q(1)) a(3)*ones(Q(1),Q(2)) a(2)*ones(Q(1),Q(3));
    a(2)*ones(Q(2),Q(1)) a(1)*ones(Q(2),Q(2)) a(3)*ones(Q(2),Q(3));
    a(3)*ones(Q(3),Q(1)) a(2)*ones(Q(3),Q(2)) a(1)*ones(Q(3),Q(3))];
%}

C_blocks=cell(M,M);
for n=1:length(Q)
    for m=1:M
        C_blocks{n,m}=c(1+mod(n-m,M))*ones(Q(n),Q(m));
    end
end
C=cell2mat(C_blocks);

W_M=(1/sqrt(M))*dftmtx(M);
W_Q=cell(1,M);
for n=1:M
    W_Q{n}=dftmtx(Q(n))*1/sqrt(Q(n)); %test sqrt and no sqrt on Q(n)
end

P_blocks=cell(M,M);
for n=1:M
    for m=1:M
        P_blocks{n,m}=W_M(n,m)*[W_Q{m}; zeros(Q(1)-Q(m),Q(m))];
    end
end

P=cell2mat(P_blocks);
pseudo_diag=P*C*P';
pseudo_diag=pseudo_diag.*(abs(pseudo_diag)>=1e-10);


Q_wave_blocks=cell(M,1);

for n=1:M
        Q_wave_blocks{n}=Q(n)*ones(Q(n),1);
end
Q_wave=cell2mat(Q_wave_blocks);


T=P*diag(a)'*diag(sqrt(Q_wave))*P';
x=P*diag(sqrt(Q_wave))*(a);
Tx=T*x;


y=ones(s,1);
non_symm_diagonal=(1/sqrt(s)).*dftmtx(s)*diag([1 1 1 1 1 1])*diag(sqrt(Q_wave))*P';
fft(Q_wave)
%[V,D]=eigs(C,[],length(Q)^2);
%D=diag(D);