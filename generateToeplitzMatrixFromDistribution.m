% Input: Q, B, rho_hat_symm_2B
% Output: The block toeplitz matrix T as shows in the second moment.
function [T] = generateToeplitzMatrixFromDistribution(Q, B, rho_hat_symm_2B)

% Assign the variations of Rho coefficients:
assignRhoCoefficients(Q, B, rho_hat_symm_2B);

% Build the Block Toeplitz matrix:
T = BTTB(repmat(rho_hat_half_2B.', Q, 1));
end