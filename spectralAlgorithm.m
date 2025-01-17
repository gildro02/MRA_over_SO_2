% This function runs the FM algorithm on the 2nd order moments in
% order to recover the approximate coefficients of the picture, and the
% corresponding error.
% Input: the moments, Q, B, sigma, the original coefficients a_symm_1B (only used
% for finding the error, not used in the algorithm itself).
% Output: The coefficients of the approximate solution (rotated by the
% correct angle) and the corresponding squared error.
% Optionally - D_sorted and kappa are also returned, where D_sorted is the
% vector of eigenvalues sorted from most positive to most negative, and
% kappa is the index of eigenvector picked by the algorithm (the
% eigenvector whos eigenvalue is most isolated).
function [a_est_rotated, error_squared_spectral, D_sorted, kappa] = spectralAlgorithm(M_1, M_2, parameters)

% Get the relevant parameters.
[Q, B, sigma, a_symm_1B, resolution_error, force_pure_phases, spectral_error_func, spectral_algorithm_version] = ...
    getStructFields(parameters, 'Q', 'B', 'sigma', 'a_symm_1B', 'resolution_error',...
    'force_pure_phases', 'spectral_error_func', 'spectral_algorithm_version');

% Assign the frequencies:
assignFrequencies(Q, B);

% Build the 2D DFT matrix:
W = DFT_Matrix_2D(Q, B);

% Debiasing
M_2_wave = M_2 - biasMatrixForSecondMoment(Q, B, sigma);
% Get the empirical power spectrum:
P_a = diag(M_2_wave);
v_p = 1 ./ sqrt(P_a);

if strcmp(spectral_algorithm_version, "old")
    % Remove power spectrum from M_2 and use inverse DFT via conjugation by W:
    M_2_mathcal = (1 / (2*pi)) * W' * diag(v_p) * M_2_wave * diag(v_p) * W;
    
    %Force Hermitian (it is slightly not hermitian because of numerics).
    M_2_mathcal = (M_2_mathcal' + M_2_mathcal) / 2;
    
    % Get the eigenvalue decomposition, sort eigenvalues from most positive to
    % most negative
    [V, D] = eig(M_2_mathcal, "vector");
    [D_sorted, I_sort] = sort(D, "descend");
    V_sorted = V(:, I_sort);
    
    % Pick eigenvalue farthest away from any other
    % eigenvalue (and its corresponding eigenvector):
    [~, kappa] = max(min(abs(D_sorted - D_sorted.') + diag(inf(1, Q * (2*B + 1))), [], 2));
    v_kappa = V_sorted(:, kappa);
    
    % Get the approximate vector of coefficient phases:
    a_est_tilde = vec(fft2(mat(v_kappa, Q, B)));
    
elseif strcmp(spectral_algorithm_version, "new")
    % Remove power spectrum from M_2 (DO NOT CONJUGATE BY W). Also, no
    % dividing by pi.
    M_2_mathcal = diag(v_p) * M_2_wave * diag(v_p);
    
    %Force Hermitian (it is slightly not hermitian because of numerics).
    M_2_mathcal = (M_2_mathcal' + M_2_mathcal)/2;
    
    % Get the eigenvalue decomposition, sort eigenvalues from most positive to
    % most negative
    [V, D] = eig(M_2_mathcal, "vector");
    [D_sorted, I_sort] = sort(D, "descend");
    V_sorted = V(:, I_sort);
    
    % Pick eigenvalue farthest away from any other
    % eigenvalue (and its corresponding eigenvector):
    [~, kappa] = max(min(abs(D_sorted - D_sorted.') + diag(inf(1, Q * (2*B + 1))), [], 2));
    
    % Get the approximate vector of coefficient phases:
    a_est_tilde = sqrt(Q * (2*B + 1)) .* V_sorted(:, kappa);
else
    error("Invalid input for spectral algorithm version");
end


if force_pure_phases == true
    a_est_tilde = a_est_tilde ./ abs(a_est_tilde);
end
% Re-introduce the power spectrum and use M_1 for phase determination via
% the 0'th frequency.
a_est = a_est_tilde...
    .* exp(1i .* (angle(M_1(k_symm_1B_Q_vec == 0 & q_full_1B_vec == 0))...
    - angle(a_est_tilde(k_symm_1B_Q_vec == 0 & q_full_1B_vec == 0))))...
    .* sqrt(P_a);

if strcmp(spectral_error_func, 'unrestricted') == 1
    % Find the error of The algorithm, and the corresponding rotation
    % angle, via the unrestricted option (default).
    [error_squared_spectral, theta_min_spectral] = ...
        circ_error_continuous_unrestricted_2D(a_est, a_symm_1B, Q, B, resolution_error);
    
elseif strcmp(spectral_error_func, 'restricted') == 1
    % Find the error of The algorithm, and the corresponding rotation
    % angle via the restricted option.
    % l/(1 + 2*B) is the fraction of rotation of the result, meaning that
    % 2*pi/(1 + 2*B) is the corresponding angle of rotation.
    [error_squared_spectral, l] = ...
        Circ_Error_Continuous_2D(a_est, a_symm_1B, Q, B);
    theta_min_spectral = 2 * pi * l / (2*B + 1);
    
else
    error("Illegal Error Function Input!")
end

% Rotate back by the rotation angle.
a_est_rotated = rotateImageViaCoefficients(Q, B, a_est, theta_min_spectral);

end