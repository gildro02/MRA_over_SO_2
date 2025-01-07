function [C] = generateCirculantMatrixFromDistribution(Q, B, rho_hat_symm_2B)

% Assign the frequencies:
assignFrequencies(Q, B);

% Assign the variations of Rho coefficients:
assignRhoCoefficients(Q, B, rho_hat_symm_2B);

% Define relevant vectors and matrices:
h = (1 / (2*B + 1)) * (k_half_2B.' .* [0; conj(flip(rho_hat_half_2B(2:end)))]...
    +(2*B + 1 - k_half_2B.') .* rho_hat_half_2B); %rho_hat_half_2B is defined in assignRhoCoefficients.
zeta = repmat(h.', [Q, 1]);
C = BCCB(zeta);

end