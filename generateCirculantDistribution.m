% Input: coefficients of a given distribution, vector t of times over which
% the values are non-negative.
% Output: distribution function and coefficient vector of distribution that
% has the circulant property.
function [rho_hat_symm_2B_circ, rho_func_circ] = generateCirculantDistribution(rho_hat_symm_2B, t, Q, B)

% Assign the frequencies:
assignFrequencies(Q, B);

% Assign the variations of Rho coefficients:
assignRhoCoefficients(Q, B, rho_hat_symm_2B);

% Find Coefficient vector that adheres to the circulant condition:
h = (1 / (2*B + 1)) * (k_half_2B.' .* [0; conj(flip(rho_hat_half_2B(2:end)))]...
    +(2*B + 1 - k_half_2B.') .* rho_hat_half_2B); %rho_hat_half_2B is defined in assignRhoCoefficients.

% Make temporary distribution, could be negative
rho_hat_symm_2B_temp = [conj(flip(h(2:end))); h];
rho_func_temp = @(x) real(rho_hat_symm_2B_temp.'*exp(1i*x'*k_symm_2B).');

% Initialize the output distribution:
rho_hat_symm_2B_circ = rho_hat_symm_2B_temp;

% Ensure the distribution is non-negative (by removing the minimum from the 0'th coefficient):
min_rho_func = min(rho_func_temp(t));
rho_hat_symm_2B_circ(k_symm_2B == 0) = rho_hat_symm_2B_temp(k_symm_2B == 0) - min(min_rho_func, 0); %make >=0

% Renormalize the coefficients s.t. the 0'th coefficient is 1/2pi:
rho_hat_symm_2B_circ = rho_hat_symm_2B_circ .* (1/(2*pi)) * (1 / rho_hat_symm_2B_circ(k_symm_2B == 0)); %renormalize

% Build the distribution function:
rho_func_circ = @(x) real(rho_hat_symm_2B_circ.'*exp(1i*x.'*k_symm_2B).'); %the final rho_func

end