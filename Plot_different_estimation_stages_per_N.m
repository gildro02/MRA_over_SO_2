%% Define variables
B = 20;
Q = 10;
t = -pi:1e-4:pi;
resolution_error = 1e-4;
moment_type = 'empirical';
%moment_type = 'analytical';
force_pure_phases = false;
spectral_algorithm_version = 'new';
sigma = 1;
% N_vec = [1, 1e3, 1e6];
N_vec = [1]
% N = 3e5;
% isUniformPowerSpectrum = true;
isUniformPowerSpectrum = false;

size_image = 101;

%% Generate distribution:
f = 0.1;
amp = 10;
% % Load temp for circulant distribution:
% load(fullfile('circulant_seeds', 'temp_circulant_random.mat'),...
%     "temp_circulant")
% % Add perturbation from circulant distribution:
% temp = temp_circulant...
%     .* exp(1i * (f * sqrt((1:2*B).'))) * amp;
temp = sqrt(1:2*B).';
[rho_hat_symm_2B, ~] = GenerateDistribution(temp, t, Q, B);
[rho_hat_symm_2B, rho_func] = generateCirculantDistribution(rho_hat_symm_2B, t, Q, B);
distanceFromCirculant(Q, B, rho_hat_symm_2B)
%% Make image:
image_folder = 'Flower_Images';
% image_name = 'Oxalis_tetraphylla_flower.jpg';
%image_name = 'flower-circle-border.jpg';
%image_name = 'nonzero_coefficients_offdiagonal.jpg';
image_name = 'fading_angular_dependence_shape.jpg';


% UNIFORM POWER SPECTRUM
[a_symm_1B, image, Phi_ns_mat] = ...
    Generate_Picture_cut(['.\' image_folder '\' image_name], B, Q, size_image,...
    isUniformPowerSpectrum);

energy_of_a = sum(abs(a_symm_1B).^2);
SNR = energy_of_a ./ ((2*B+1) * Q * sigma .^ 2);
%% define params vector
diff_N_vec = [N_vec(1) diff(N_vec)];
for n_idx = 1:length(N_vec)
    params_vec_diff_N(n_idx).N = diff_N_vec(n_idx);
    params_vec_diff_N(n_idx).Q = Q;
    params_vec_diff_N(n_idx).B = B;
    params_vec_diff_N(n_idx).sigma = sigma;
    params_vec_diff_N(n_idx).a_symm_1B = a_symm_1B;
    params_vec_diff_N(n_idx).rho_func = rho_func;
    params_vec_diff_N(n_idx).t = t;
    params_vec_diff_N(n_idx).resolution_error = resolution_error;
    % a_symm_1B_rot = GetCoeffRotated(params_for_single_picture);
end

%% Generate Moments & make images:
a_approx_symm_1B_running = zeros(length(a_symm_1B), length(N_vec));
error_squared_FM = zeros(length(N_vec), 1);
M_1_running = zeros(length(a_symm_1B), length(N_vec));
M_2_running = zeros(length(a_symm_1B), length(a_symm_1B), length(N_vec));
for n_idx = 1:length(N_vec)
    [M_1_partial, M_2_partial] = generateEmpiricalMoments(params_vec_diff_N(n_idx));
        if n_idx == 1
            M_1_running(:, n_idx) = M_1_partial;
            M_2_running(:, :, n_idx) = M_2_partial;
        else
            M_1_running(:, n_idx) = (N_vec(n_idx - 1) .* M_1_running(:, n_idx - 1) + diff_N_vec(n_idx) .* M_1_partial) ./ N_vec(n_idx);
            M_2_running(:, :, n_idx) = (N_vec(n_idx - 1) .* M_2_running(:, :, n_idx - 1) + diff_N_vec(n_idx) .* M_2_partial) ./ N_vec(n_idx);
        end
%     M_1_running(:, n_idx) = M_1_partial;
%     M_2_running(:, :, n_idx) = M_2_partial;
    [a_approx_symm_1B_running(:, n_idx), error_squared_FM(n_idx)] =...
        frequencyMarchingAlgorithm(M_1_running(:, n_idx), M_2_running(:, :, n_idx), params_vec_diff_N(n_idx));
end

%% Plot the Pictures
close all
fig_true = PlotImageOfCoeff(a_symm_1B, Phi_ns_mat, size_image, false);
for n_idx = 1:length(N_vec)
    fig_partial(n_idx) = PlotImageOfCoeff(a_approx_symm_1B_running(:, n_idx), Phi_ns_mat, size_image, false);
    
end

%% Save plots
save_path = fullfile(".", "Figures_Thesis", "Example_Figs_Convergence");
load(fullfile(save_path, "attempt_idx.mat"));
attempt_idx = attempt_idx + 1;
save_path_attempt = fullfile(save_path, "attempt_" + attempt_idx);
mkdir(save_path_attempt);

saveFigureToAllFormats(fig_true, "fig_true", save_path_attempt);
for n_idx = 1:length(N_vec)
    saveFigureToAllFormats(fig_partial(n_idx), "fig_N_is_" + N_vec(n_idx), save_path_attempt);
end

sizeThresholdKB = 2000; %CHECK BEFORE LARGE RUN!!!!
saveNumericVariablesBelowThreshold(save_path_attempt, "all_numeric_variables", sizeThresholdKB);

save(fullfile(save_path, "attempt_idx.mat"), 'attempt_idx')
