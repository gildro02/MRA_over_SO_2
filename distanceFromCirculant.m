% Input: Q, B, and the distribution coefficients
% Output: The squared distance of T from its circulant approximation C, equivalent to
% the distance between M_2_C_mathcal and M_2_T_mathcal.
function dist_from_circ = distanceFromCirculant(Q, B, rho_hat_symm_2B)

% Assign the Rho coefficients and frequencies:
assignFrequencies(Q, B);
assignRhoCoefficients(Q, B, rho_hat_symm_2B);

% The formula squared distance of T from its circulant approximation C. This is the exact
% formula for the distance, as written in the paper:
dist_from_circ = (Q ^ 2) * sum(abs((conj(rho_hat_half_2B(2:end)) - flip(rho_hat_half_2B(2:end))) .^ 2)...
    .*(k_half_2B(2:end).' .* (flip(k_half_2B(2:end).'))) ./ (2*B + 1));
end