% Input: Q, B, sigma
% Output: the bias matrix for this noise level.
function [bias_mat] = biasMatrixForSecondMoment(Q, B, sigma)

% Assign the frequencies:
assignFrequencies(Q, B);

% Bias occurs only for q1 = q2, k1 = +-k2
Logical_Bias_Mat = zeros((2*B + 1) * Q);
Logical_Bias_Mat(k_symm_1B_Q_vec - k_symm_1B_Q_vec.' == 0&...
    q_full_1B_vec - q_full_1B_vec.' == 0) = 1; %q1 = q2, k1 = k2
Logical_Bias_Mat(k_symm_1B_Q_vec + k_symm_1B_Q_vec.' == 0&...
    q_full_1B_vec - q_full_1B_vec.' == 0) = 1; %q1 = q2, k1 = -k2

bias_mat = (sigma.^2) * Logical_Bias_Mat;
end