function [bound] = boundForSpectralAlgorithm_NewVersion(varargin)

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

% Run the spectral algorithm to get the index of the eigenvector picked by 
% the algorithm:
[~, ~, ~, kappa_T] = spectralAlgorithm(M_1_T, M_2_T, parameters);

% Generate T and C to find their eigenvalues:
T = generateToeplitzMatrixFromDistribution(Q, B, rho_hat_symm_2B);
C = generateCirculantMatrixFromDistribution(Q, B, rho_hat_symm_2B);

[~, D_T] = eig(T, "vector");
D_T_sorted = sort(D_T, "descend");
[~, D_C] = eig(C, "vector");
D_C_sorted = sort(D_C, "descend");


% Get the minimum distance between the picked true eigenvalue and the
% corresponding eigenvalue of the circulant case (or, the other way around
% - the maximum between them):
delta_kappa = max(min(abs(D_T_sorted(kappa_T) - D_C_sorted(1:end ~= kappa_T))), ...
    min(abs(D_C_sorted(kappa_T) - D_T_sorted(1:end ~= kappa_T))));

% Get the squared distance of T from its circulant approximation C, equivalent to
% the distance between M_2_C_mathcal and M_2_T_mathcal.
dist_from_circ = distanceFromCirculant(Q, B, rho_hat_symm_2B);

% Get the bound on the squared error of the spectral algorithm:
bound = 2 * Q * (2*B + 1) * max(P_a) .* (1 - sqrt(1 - dist_from_circ ./ delta_kappa .^ 2));
end