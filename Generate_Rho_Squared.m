function [rho_squared_coeff,rho_squared]=Generate_Rho_Squared(seed_rho_coeff)
%returns positive, normilized, Band-Limited distribution from seed
%coefficients seed_rho_coeff. Distribution is periodic on [0,2pi];

B=length(seed_rho_coeff);
s=size(seed_rho_coeff);
if s(2)==1
    seed_rho_coeff=seed_rho_coeff.';
end
rho_coeff=[conj(flip(seed_rho_coeff)) 1/(2*pi) seed_rho_coeff]; %1/2pi for normilization of distribution to 1

%{
full_rho_freq=-B:B;
rho=@(x) real(rho_coeff*exp(1i*x'*full_rho_freq)');
%}
%{
Norm_Rho_Squared=integral(@(x) rho(x).^2,0,2*pi);
rho_squared_analitic=@(x) (rho(x).^2); %./Norm_Rho_Squared
%}
rho_squared_coeff=conv(rho_coeff,rho_coeff);%*(1/(2*pi)); %ensure positivity
rho_squared_coeff=rho_squared_coeff./(2*pi*rho_squared_coeff(2*B+1)); %normilize
full_rho_squared_freq=-2*B:2*B;
rho_squared=@(x) real(rho_squared_coeff*exp(1i*x'*full_rho_squared_freq).');

if s(2)==1
    rho_squared_coeff=rho_squared_coeff.';
end

end