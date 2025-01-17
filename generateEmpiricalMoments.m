% Input: a struct of all relavent parameters (N, B, Q, t, rho_func, sigma,
% a_symm_1B)
% Output: the first and second empirical moments.
function [M_1_emp, M_2_emp] = generateEmpiricalMoments(parameters)

% Add all the relevent parameters to the workspace:
[N, Q, B, sigma, a_symm_1B, rho_func, t] = getStructFields(parameters, 'N', 'Q', 'B', 'sigma', 'a_symm_1B', 'rho_func', 't');

% Assign all the frequency vectors:
assignFrequencies(Q,B);

% Assign all the relavent orientations of a_symm_1B:
a_symm_1B_mat = mat(a_symm_1B, Q, B);
a_half_1B_mat = a_symm_1B_mat(:, k_symm_1B >= 0);
a_half_1B_vec = vec(a_half_1B_mat);

% Randomize N rotations from distribution rho_func, over possible
% rotations t:
rotations = randsample(t, N, true, abs(rho_func(t))).';

% Generate rotated versions of a_half_1B_mat:
% real noise:
%a_half_1B_mat_rot = a_half_1B_mat .* exp(-1i * k_half_1B_Q_mat .* permute(rotations, [3, 2, 1]))...
%    + sigma * randn([size(a_half_1B_mat) N]);

% Complex noise:

% Noise at non-zero frequencies - var = sigma^2/2 for both real and
% imaginary part.
size_a_half_1B_mat_positive_k = [Q, B];
noise_positive_freq = (sigma / sqrt(2)) .* (randn([size_a_half_1B_mat_positive_k N])...
    + 1i .* randn([size_a_half_1B_mat_positive_k N]));

% Noise at frequency zero - var = sigma^2, pure real number.
noise_0_freq = sigma .* randn([Q 1 N]);

% Combined noise matrix:
noise = cat(2, noise_0_freq, noise_positive_freq);

% Rotate Images and add noise:
a_half_1B_mat_rot = a_half_1B_mat .* exp(-1i * k_half_1B_Q_mat .* permute(rotations, [3, 2, 1]))...
    + noise;

% Define the full (rotated versions of) a_symm_1B_rot via conjucation and
% flip:
a_half_1B_mat_rot_conj = conj(flip(reshape(...
    a_half_1B_mat_rot(repmat(k_half_1B_Q_mat ~= 0, [1, 1, N]))...
    , [Q, B, N]), 2));
a_symm_1B_mat_rot = cat(2, a_half_1B_mat_rot_conj, a_half_1B_mat_rot);
a_symm_1B_rot = reshape(a_symm_1B_mat_rot, [Q * (2*B+1), N]);

% Generate empirical moments:
M_1_emp = vec(mean(a_symm_1B_mat_rot, 3));
M_2_emp = (1/N) * (a_symm_1B_rot * a_symm_1B_rot'); %'=complex conjugate transpose
end