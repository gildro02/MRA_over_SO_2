%% Define variables
N = 1e6;
B = 10;
Q = 2;
t = -pi:1e-4:pi;
resolution_error = 1e-4;
moment_type = 'empirical';
%moment_type = 'analytical';
force_pure_phases = false;

num_unique_sigma = 80;
num_rep_sigma = 400;
sigma_vec_reduced = logspace(-1.5, 0.5, num_unique_sigma).';
%sigma_vec_reduced = logspace(-1.5, 1.2, num_unique_sigma).';
sigma_vec = repelem(sigma_vec_reduced, num_rep_sigma);

isUniformPowerSpectrum = true;
%% Generate distribution:
f = 0.1;
amp = 10;
% Load temp for circulant distribution:
load(fullfile('circulant_seeds', 'temp_circulant_random.mat'),...
    "temp_circulant")
% Add perturbation from circulant distribution:
temp = temp_circulant...
    .* exp(1i * (f * sqrt((1:2*B).'))) * amp;

% Generate Distribution:
[rho_hat_symm_2B, rho_func] = GenerateDistribution(temp, t, Q, B);
dist_from_circ = distanceFromCirculant(Q, B, rho_hat_symm_2B);

%% Make image:
size_image = 51;
image_folder = 'Flower_Images';
image_name = 'Oxalis_tetraphylla_flower.jpg';
[a_symm_1B, image, Phi_ns_mat] = ...
    Generate_Picture_cut(['.\' image_folder '\' image_name], B, Q, size_image,...
    isUniformPowerSpectrum);

%% Initialize error vectors:
error_squared_FM = zeros(length(sigma_vec), 1);
error_squared_spectral = zeros(length(sigma_vec), 1);
M_1_error = zeros(length(sigma_vec), 1);
M_2_error = zeros(length(sigma_vec), 1);
%% Run the experiment
tic
for n = 1:length(sigma_vec)
    [error_squared_FM(n), coeff_FM, error_squared_spectral(n), coeff_spectral, M_1_error(n), M_2_error(n)] = ...
        runNumericalExperiment('B', B, 'Q', Q, 'N', N, 'sigma', sigma_vec(n), 'rho_func',...
        rho_func, 't', t, 'rho_hat_symm_2B', rho_hat_symm_2B,...
        'a_symm_1B', a_symm_1B, 'moment_type', moment_type,...
        'resolution_error', resolution_error, 'force_pure_phases', force_pure_phases);
    disp("n = " + n);
    toc
end

%% Plot the errors

energy_of_a = sum(abs(a_symm_1B) .^ 2);
SNR = energy_of_a ./ ((2*B + 1) * Q * sigma_vec .^ 2);

error_squared_spectral_relative = error_squared_spectral ./ energy_of_a;
error_squared_FM_relative = error_squared_FM ./ energy_of_a;
prctile_region = 20;
asymp_start_percent = 50;

% Get asymptote percent:
min_err_spectral = ...
    asymptoticSpectralError(Q, B, a_symm_1B, rho_hat_symm_2B, resolution_error);
min_err_spectral_relative = min_err_spectral ./ energy_of_a;

fig = plotMRAError('variable', SNR, 'num_rep', num_rep_sigma, 'mean_or_median', 'median',...
    'error_spectral', error_squared_spectral_relative, 'error_FM', error_squared_FM_relative,...
    'labels', {"SNR", "Relative Squared Error"}, 'title', '',...
    'std_factor_or_prctile_region', prctile_region, 'is_x_math', false,...
    'y_asymptote_start_percent', asymp_start_percent, 'y_asymptote_val', min_err_spectral_relative);

%% Save the plot:

% Push test_number up by 1, make directory of test
dir_path = makeNewTestDir();
figure_name = "func_of_sigma_Uniform_PS";
saveFigureToAllFormats(fig, figure_name, dir_path);
sizeThresholdKB = 2000;
saveNumericVariablesBelowThreshold(dir_path, "all_numeric_variables", sizeThresholdKB);
