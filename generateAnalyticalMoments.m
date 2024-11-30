% Input: a struct of all relavent parameters (B, Q, rho_hat_symm_2B, sigma,
% a_symm_1B)
% Output: the first and second analytical moments.
function [M_1_true, M_2_true] = generateAnalyticalMoments(parameters)

% Add all the relevent parameters to the workspace:
[B, Q, sigma, a_symm_1B, rho_hat_symm_2B] = getStructFields(parameters, 'B', 'Q', 'sigma', 'a_symm_1B', 'rho_hat_symm_2B');

% Assign all the frequency vectors, and all the rho coefficient vectors:
assignFrequencies(Q, B);
assignRhoCoefficients(Q, B, rho_hat_symm_2B);

%% 1st Moment:
M_1_true = a_symm_1B .* varrho .* (2*pi);

%% 2nd Moment:
T = BTTB(repmat(rho_hat_half_2B.', Q, 1));

% Bias Matrix:
[bias_mat] = biasMatrixForSecondMoment(Q, B, sigma);

M_2_true =(2*pi) * diag(a_symm_1B) * T * diag(conj(a_symm_1B)) + bias_mat; %'=complex conjugate transpose

end