% This function runs the FM algorithm on the 2nd order moments in
% order to recover the approximate coefficients of the picture, and the
% corresponding error.
% Input: the moments, Q, B, sigma, the original coefficients a_symm_1B (only used
% for finding the error, not used in the algorithm itself).
% resolution_error_FM is the resolution for determination of the FM error
% (the size of d_theta that partitions [0, 2pi]).
% Output: The coefficients of the approximate solution (rotated by the
% correct angle) and the corresponding squared error.
function [a_approx_symm_1B_rotated, error_squared_FM] = frequencyMarchingAlgorithm(M_1, M_2, parameters)

[Q, B, sigma, a_symm_1B, resolution_error] = getStructFields(parameters, 'Q', 'B', 'sigma', 'a_symm_1B', 'resolution_error');

% Assign the frequencies:
assignFrequencies(Q, B);

% Debiasing:
M_2_wave = M_2 - biasMatrixForSecondMoment(Q, B, sigma);
S_mathcal = 2 * pi * diag(1./M_1) * M_2_wave * diag(1./M_1)';

%Build weights per each block of constant k:
weights_q = cell(2*B + 1, 2*B + 1);

for k1 = -B:B
    for k2 = -B:B
        weight = exp(-q_1B).' * exp(-q_1B); %low q's are worth most
        weights_q{k_symm_1B == k1, k_symm_1B == k2} = weight ./ sum(weight, "all"); %normalize
    end
end

S_mathcal_tilde = zeros(2*B + 1, 2*B + 1); %condensed S
for k1 = -B:B
    for k2 = -B:B
        S_mathcal_tilde(k_symm_1B == k1, k_symm_1B == k2)...
            = sum(weights_q{k_symm_1B == k1, k_symm_1B == k2}...
            .* S_mathcal(k_symm_1B_Q_vec == k1, k_symm_1B_Q_vec == k2), "all");
    end
end

weights_k = cell(1, B);
for k = 2:B
    weights_k{k} = exp(-abs(ceil(k/2) - (1:k-1).')); %lorentzian of length k-1 centered at ceil(length/2)
    weights_k{k} = weights_k{k} ./ sum(weights_k{k}); %normalize to sum 1.
end

rho_approx_half_1B = zeros(length(k_half_1B), 1);
rho_approx_tilde_half_1B = zeros(length(k_half_1B), 1);

rho_abs_half_1B = sqrt(abs(1 ./ (2 * pi * diag(S_mathcal_tilde(k_symm_1B >= 0, k_symm_1B >= 0)))));
rho_approx_half_1B(k_half_1B == 0) = 1/(2*pi);
%arbitrarily pick phase of 1'st coeff to be 1:
rho_approx_half_1B(k_half_1B == 1) = rho_abs_half_1B(k_half_1B == 1);

for k=2:B
    rho_approx_tilde_half_1B(k_half_1B == k) = sum(weights_k{k} .* flip(rho_approx_half_1B(1 <= k_half_1B & k_half_1B <= k-1))./...
        conj(rho_approx_half_1B(1 <= k_half_1B & k_half_1B <= k-1)) ./ ...
        S_mathcal_tilde(k_symm_1B == k, 1 <= k_symm_1B & k_symm_1B <= k-1).');
    
    rho_approx_half_1B(k_half_1B == k) = rho_abs_half_1B(k_half_1B == k) .* ...
        rho_approx_tilde_half_1B(k_half_1B == k) ./ abs(rho_approx_tilde_half_1B(k_half_1B == k));
end
rho_approx_symm_1B = [conj(flip(rho_approx_half_1B(2:end))); rho_approx_half_1B];
varrho_approx = repelem(rho_approx_symm_1B, Q);

a_approx_symm_1B = M_1 ./ (2 * pi .* varrho_approx);

[error_squared_FM, theta_rotation]=...
    circ_error_continuous_unrestricted_2D(a_approx_symm_1B, a_symm_1B, Q, B, resolution_error);

% theta_rotation is the relative rotation between the approximate picture
% and the true picture.
a_approx_symm_1B_rotated = rotateImageViaCoefficients(Q, B, a_approx_symm_1B, theta_rotation);
end