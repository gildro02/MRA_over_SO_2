% Input: a struct of all relavent parameters (B, Q, rho_hat_symm_2B, sigma,
% a_symm_1B)
% Output: the second Biased Circulant moment.
function [M_2_C] = generateCirculantSecondMomentBiased(parameters)

% Add all the relevent parameters to the workspace:
[Q, B, sigma, a_symm_1B, rho_hat_symm_2B] = ...
    getStructFields(parameters, 'Q', 'B', 'sigma', 'a_symm_1B', 'rho_hat_symm_2B');

% Assign the frequencies:
assignFrequencies(Q, B);

% Assign the variations of Rho coefficients:
assignRhoCoefficients(Q, B, rho_hat_symm_2B);

% Define relevant vectors and matrices:
h = (1 / (2*B + 1)) * (k_half_2B.' .* [0; conj(flip(rho_hat_half_2B(2:end)))]...
    +(2*B + 1 - k_half_2B.') .* rho_hat_half_2B); %rho_hat_half_2B is defined in assignRhoCoefficients.
zeta = repmat(h.', [Q, 1]);
C = BCCB(zeta);

% Define the unbiased circulant moment + Add the bias:
M_2_C = 2 * pi * diag(a_symm_1B) * C * diag(a_symm_1B)'...
    + biasMatrixForSecondMoment(Q, B, sigma);

end