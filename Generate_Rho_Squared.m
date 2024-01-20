% Input: 1) seed: seed coefficients of the resulting distribution,
% corresponding to frequencies 1:B (a vector of length 1B).
% Output: 1) rho_squared_coeff: the 4B+1 x 1 vector of Fourier coefficients of the
% distribution rho_squared, corresponding to frequencies -2B:2B.
% 2) rho_squared: the resulting distribution, defined by first completeing
% 'seed' to frequencies -B:B (with 1/2pi as a middle value), then squaring
% the function by convoluting this vector with itself (this ensures
% non-negativity), and then renormalizing.

function [rho_squared_coeff,rho_squared]=Generate_Rho_Squared(seed)
%returns positive, normilized, Band-Limited distribution from seed
%coefficients seed. Distribution is periodic on [0,2pi];

B=length(seed);
s=size(seed);
if s(2)==1
    seed=seed.';
end
rho_coeff=[conj(flip(seed)) 1/(2*pi) seed]; %1/2pi for normilization of distribution to 1

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