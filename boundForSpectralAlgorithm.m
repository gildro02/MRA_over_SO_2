function [bound] = boundForSpectralAlgorithm(varargin)

% Create an input parser object
p = generateMRAInputParser();

% Parse the paramters
parse(p, varargin{:})
parameters = p.Results;

[Q, B, a_symm_1B, rho_hat_symm_2B] = ...
    getStructFields(parameters, 'Q', 'B', 'a_symm_1B', 'rho_hat_symm_2B');

% Assign all the frequency vectors, and all the rho coefficient vectors:
assignFrequencies(Q, B);
assignRhoCoefficients(Q, B, rho_hat_symm_2B);

% Get the power Spectrum of a_symm_1B:
P_a = abs(a_symm_1B) .^ 2;

% Define both the true moment and the moment with the approximated
% circulant matrix:
M_2_C = generateCirculantSecondMomentBiased(parameters);
[M_1_T, M_2_T] = generateAnalyticalMoments(parameters);

% Run the spectral algorithm on both, extracting the eigenvalues and index
% of picked eigenvector:
[~, ~, D_T_sorted, kappa_T] = spectralAlgorithm(M_1_T, M_2_T, parameters);
[~, ~, D_C_sorted, ~] = spectralAlgorithm(M_1_T, M_2_C, parameters);

% Get the minimum distance between the picked true eigenvalue and the
% corresponding eigenvalue of the circulant case (or, the other way around
% - the maximum between them):
delta_kappa = max(min(abs(D_T_sorted(kappa_T) - D_C_sorted(1:end ~= kappa_T))), ...
    min(abs(D_C_sorted(kappa_T) - D_T_sorted(1:end ~= kappa_T))));

% Get the distance of T from its circulant approximation C, equivalent to
% the distance between M_2_C_mathcal and M_2_T_mathcal. This is the exact
% formula for the distance, as written in the paper.
dist_from_circ = (Q ^ 2) * sum(abs((conj(rho_hat_half_2B(2:end)) - flip(rho_hat_half_2B(2:end))) .^ 2)...
    .*(k_half_2B(2:end).' .* (flip(k_half_2B(2:end).'))) ./ (2*B + 1));

% Get the bound on the squared error of the spectral algorithm:
bound = 2 * Q * (2*B + 1) * max(P_a) .* (1 - sqrt(1 - dist_from_circ ./ delta_kappa .^ 2));
end