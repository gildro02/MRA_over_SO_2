% Input: Image coefficients and distribution coefficients.
% Output: The theoretical bound for the spectral algorithm. Takes the
% minimum over rotations of rho if 'bound_rotation_resolution' is
% specified, otherwise returns the default bound.
function [bound, theta_opt, dist_from_circ_opt, delta_kappa_opt] = boundForSpectralAlgorithm(varargin)

% Create an input parser object
p = generateMRAInputParser();

% Parse the paramters
parse(p, varargin{:})
parameters = p.Results;

[Q, B, a_symm_1B, rho_hat_symm_2B, bound_rotation_resolution] = ...
    getStructFields(parameters, 'Q', 'B', 'a_symm_1B', 'rho_hat_symm_2B', 'bound_rotation_resolution');

% Assign all the frequency vectors, and all the rho coefficient vectors:
assignFrequencies(Q, B);
assignRhoCoefficients(Q, B, rho_hat_symm_2B);

% Get the power Spectrum of a_symm_1B:
P_a = abs(a_symm_1B) .^ 2;

% Define both the true moment and the moment with the approximated
% circulant matrix:
[M_1_T, M_2_T] = generateAnalyticalMoments(parameters);

% Run the spectral algorithm to get the index of the eigenvector picked by
% the algorithm:
[~, ~, ~, kappa_T] = spectralAlgorithm(M_1_T, M_2_T, parameters);


% Initialize vectors of rotations and resulting bounds:
if isempty(bound_rotation_resolution)
    theta_rot = 0;
else
    theta_rot = (0:bound_rotation_resolution:2*pi - bound_rotation_resolution).';
end
bounds = zeros(length(theta_rot), 1);
dist_from_circ = zeros(length(theta_rot), 1);
delta_kappa = zeros(length(theta_rot), 1);

% Find the bound for each theta_rot:
for theta_idx = 1:length(theta_rot)
    
    exp_multi_matrix = exp(1i .* k_symm_2B.' .* theta_rot(theta_idx));
    rho_hat_symm_2B_rot = exp_multi_matrix .* rho_hat_symm_2B;
    
    % Generate T and C to find their eigenvalues:
    T = generateToeplitzMatrixFromDistribution(Q, B, rho_hat_symm_2B_rot);
    C = generateCirculantMatrixFromDistribution(Q, B, rho_hat_symm_2B_rot);
    
    [~, D_T] = eig(T, "vector");
    D_T_sorted = sort(D_T, "descend");
    [~, D_C] = eig(C, "vector");
    D_C_sorted = sort(D_C, "descend");
    
    
    % Get the minimum distance between the picked true eigenvalue and the
    % corresponding eigenvalue of the circulant case (or, the other way around
    % - the maximum between them):
    delta_kappa(theta_idx) = max(min(abs(D_T_sorted(kappa_T) - D_C_sorted(1:end ~= kappa_T))), ...
        min(abs(D_C_sorted(kappa_T) - D_T_sorted(1:end ~= kappa_T))));
    
    % Get the squared distance of T from its circulant approximation C, equivalent to
    % the distance between M_2_C_mathcal and M_2_T_mathcal.
    dist_from_circ(theta_idx) = distanceFromCirculant(Q, B, rho_hat_symm_2B_rot);
    
    % Get the bound on the squared error of the spectral algorithm:
    bounds(theta_idx) = 2 * Q * (2*B + 1) * max(P_a) .* (1 - sqrt(...
        1 - dist_from_circ(theta_idx) ./ delta_kappa(theta_idx) .^ 2));
end

% Find the minimal bound. If the resolution was not specified, "bounds" has
% 1 element and theta_min = 0;
[bound, idx_min] = min(bounds);
theta_opt = theta_rot(idx_min);
dist_from_circ_opt = dist_from_circ(idx_min);
delta_kappa_opt = delta_kappa(idx_min);
end