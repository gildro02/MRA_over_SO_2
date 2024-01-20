function [rho_hat_symm_2B, rho_func]=GenerateDistribution(seed,t,Q,B)
% Generates a legitemate distribution from a given coefficient vector.

[freq_file_name]=saveFrequencies(Q,B);
load(freq_file_name);

rho_hat_symm_2B_non_normalized =[conj(flip(seed)); 1/(2*pi); seed];
rho_hat_symm_2B_temp=(1/(2*pi)) * (1/abs(rho_hat_symm_2B_non_normalized(k_symm_2B==0)))...
    * rho_hat_symm_2B_non_normalized;
rho_func_temp = @(x) real(rho_hat_symm_2B_temp.'*exp(1i*x'*k_symm_2B).');

min_rho_func=min(rho_func_temp(t));
rho_hat_symm_2B=rho_hat_symm_2B_temp;
rho_hat_symm_2B(k_symm_2B==0)=rho_hat_symm_2B(k_symm_2B==0)-min(min_rho_func,0); %make >=0
rho_hat_symm_2B=rho_hat_symm_2B.*(1/(2*pi))*(1/rho_hat_symm_2B(k_symm_2B==0)); %renormalize
rho_func = @(x) real(rho_hat_symm_2B.'*exp(1i*x.'*k_symm_2B).'); %the final rho_func

end