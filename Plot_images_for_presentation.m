%% Define variables
B = 20;
Q = 10;
t = -pi:1e-4:pi;
resolution_error = 1e-4;
moment_type = 'empirical';
%moment_type = 'analytical';
force_pure_phases = false;
spectral_algorithm_version = 'new';
sigma = 0.2;
N = 3e5;
% N=1e5;
% N = 1e6;
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
% image_name = 'fading_angular_dependence_shape.jpg';
% image_name = 'strong_asymmetric_shape.jpg';
image_name = 'yin_yang_style.png';

% UNIFORM POWER SPECTRUM
[a_symm_1B, image, Phi_ns_mat] = ...
    Generate_Picture_cut(['.\' image_folder '\' image_name], B, Q, size_image,...
    isUniformPowerSpectrum);

energy_of_a = sum(abs(a_symm_1B).^2);
SNR = energy_of_a ./ ((2*B+1) * Q * sigma .^ 2);
%% define params vector
params_for_single_picture.N = 1;
params_for_single_picture.Q = Q;
params_for_single_picture.B = B;
params_for_single_picture.sigma = sigma;
params_for_single_picture.a_symm_1B = a_symm_1B;
params_for_single_picture.rho_func = rho_func;
params_for_single_picture.t = t;
a_symm_1B_rot = GetCoeffRotated(params_for_single_picture);
%% Run the experiment
tic
[error_squared_FM, coeff_FM, error_squared_spectral, coeff_spectral, M_1_error, M_2_error] = ...
    runNumericalExperiment('B', B, 'Q', Q, 'N', N, 'sigma', sigma, 'rho_func',...
    rho_func, 't', t, 'rho_hat_symm_2B', rho_hat_symm_2B,...
    'a_symm_1B', a_symm_1B, 'moment_type', moment_type,...
    'resolution_error', resolution_error, 'force_pure_phases', force_pure_phases,...
    'spectral_algorithm_version', spectral_algorithm_version);
toc
run_time = toc;
%% Plot the Pictures
close all


fig_true = PlotImageOfCoeff(a_symm_1B, Phi_ns_mat, size_image, false);
% title("true")
fig_noisy = PlotImageOfCoeff(a_symm_1B_rot, Phi_ns_mat, size_image, true);
% title("noisy")
fig_recovered_FM = PlotImageOfCoeff(coeff_FM, Phi_ns_mat, size_image, false);
% title("FM recovery")
fig_recovered_Spectral = PlotImageOfCoeff(coeff_spectral, Phi_ns_mat, size_image, false);
% title("spectral recovery")

%% Save plots
save_path = fullfile(".", "Figures_Thesis", "Example_Figs_For_Presentation");
load(fullfile(save_path, "attempt_idx.mat"));
attempt_idx = attempt_idx + 1;
save_path_attempt = fullfile(save_path, "attempt_" + attempt_idx);
mkdir(save_path_attempt);
saveFigureToAllFormats(fig_true, "fig_true", save_path_attempt);
saveFigureToAllFormats(fig_noisy, "fig_noisy", save_path_attempt);
saveFigureToAllFormats(fig_recovered_FM, "fig_recovered_FM", save_path_attempt);
saveFigureToAllFormats(fig_recovered_Spectral, "fig_recovered_Spectral", save_path_attempt);


sizeThresholdKB = 2000; %CHECK BEFORE LARGE RUN!!!!
saveNumericVariablesBelowThreshold(save_path_attempt, "all_numeric_variables", sizeThresholdKB);



save(fullfile(save_path, "attempt_idx.mat"), 'attempt_idx')
