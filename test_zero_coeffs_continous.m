B=2;
omega=(2*B+1);
%omega=B+1;
a_plus=1:B;
a_0=0.5;
a=[conj(flip(a_plus)) a_0 a_plus];
k=-B:B;
l=0:omega-1;
f=@(x) a*exp(1i*k.'.*x);
x=linspace(0,2*pi,omega);
f_x=f(x);
fft_f_x=fft(f_x)
sum(abs(fft_f_x)<1e-10)

freq=k.*l.';
A=exp(1i*freq*(2*pi/omega));
mat=dftmtx(omega)*A;
mat=mat.*(abs(mat)>1e-10)

dftmtx(omega)*(A*a.')