%{
v=[1 1i*rand(1,100)];
A=toeplitz(v);
B=toeplitz([0 conj(flip(v(2:end)))]);
C=[A B
    B A];
iscirculant(C)
ishermitian(C)
%}
close all
n=2000;
k=10;
r=k:-1:1;
theta=1i*2*pi./(k:-1:1);
v=r.*exp(theta);
%v=[2 1];
v(1)=abs(v(1));
T=sptoeplitz([v zeros(1,n-length(v))]);
C=sptoeplitz([v zeros(1,n-2*length(v)+1) conj(flip(v(2:end)))]);
[V_T,e_T]=eigs(T,n);
[V_C,e_C]=eigs(C,n);
eig_T=diag(e_T);
eig_C=diag(e_C);

figure
histogram(eig_T,100,"normalization","pdf")
hold on
histogram(eig_C,100,"normalization","pdf")
W=(1/sqrt(n))*dftmtx(n);
norm(W*T*W'-W*C*W',"fro")

[sorted_eig_C,indecies]=sort(eig_C);
sorted_V_C=V_C(:,indecies);
sorted_V_T=V_T(:,indecies);
abs_sorted_V_C=abs(sorted_V_C);
abs_sorted_V_T=abs(sorted_V_T);
prod=sorted_V_C./sorted_V_T;
%{
histogram(real(fft(full(T(:,1)))),100,"normalization","pdf")
histogram(abs(fft(full(T(:,1)))),100,"normalization","pdf")
histogram(imag(fft(full(T(:,1)))),100,"normalization","pdf")
legend("true eig T","true eig C","real fft T","abs fft T","imag fft T");
iscirculant(C)
%}