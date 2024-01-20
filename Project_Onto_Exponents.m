function [x_proj,coeff]=Project_Onto_Exponents(x,B)
%takes a function x and bandwidth B, returns a function x_proj that is x
%projected onto exp(-iBt),...,exp(iBt) (2B+1 functions). also returns the
%non-negative B+1 cofficients.

coeff=zeros(B+1,1);
for n=0:B
    coeff(n+1)=(1/pi^2)*(1/2*pi)*integral(@(t) x(t).*exp(-1i.*n.*t),0,2*pi);
end

x_proj=@(t) [flip(coeff(2:end))' coeff.']*exp(1i*(-B:1:B).'*t);  %also works for row vectors t, not colomns.