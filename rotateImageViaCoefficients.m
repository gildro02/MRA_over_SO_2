% Rotates a picture by an angle.
% Input: the coefficients of the picture to rotate, the angle of rotation, Q and B.
% Output: the coefficients of the rotated picture.
function [a_symm_1B_rotated] = rotateImageViaCoefficients(Q, B, a_symm_1B, theta_rotation)
phase_matrix = repmat(exp(1i .* (-B:B) .* theta_rotation), [Q, 1]);
a_symm_1B_rotated = vec(phase_matrix .* mat(a_symm_1B, Q, B));
end