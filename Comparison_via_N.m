%% Define variables
B = 10;
Q = 2;
t = -pi:1e-4:pi;
resolution_error = 1e-4;
moment_type = 'empirical';
%moment_type = 'analytical';
force_pure_phases = false;
spectral_algorithm_version = 'new';
sigma = 0.1;
num_unique_N = 120;
num_rep_N = 800;
N_vec_reduced = ceil(logspace(3, 6, num_unique_N).');
N_vec = repelem(N_vec_reduced, num_rep_N);

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

[rho_hat_symm_2B, rho_func] = GenerateDistribution(temp, t, Q, B);
dist_from_circ = distanceFromCirculant(Q, B, rho_hat_symm_2B);

%% Make image:
size_image = 51;
image_folder = 'Flower_Images';
image_name = 'Oxalis_tetraphylla_flower.jpg';
%image_name = 'flower-circle-border.jpg';
%image_name = 'nonzero_coefficients_offdiagonal.jpg';

% UNIFORM POWER SPECTRUM
[a_symm_1B, image, Phi_ns_mat] = ...
    Generate_Picture_cut(['.\' image_folder '\' image_name], B, Q, size_image,...
    isUniformPowerSpectrum);
%% Initialize error vectors:
error_squared_FM = zeros(length(N_vec), 1);
error_squared_spectral = zeros(length(N_vec), 1);
M_1_error = zeros(length(N_vec), 1);
M_2_error = zeros(length(N_vec), 1);
%% Run the experiment
tic
for n = 1:length(N_vec)
    [error_squared_FM(n), coeff_FM, error_squared_spectral(n), coeff_spectral, M_1_error(n), M_2_error(n)] = ...
        runNumericalExperiment('B', B, 'Q', Q, 'N', N_vec(n), 'sigma', sigma, 'rho_func',...
        rho_func, 't', t, 'rho_hat_symm_2B', rho_hat_symm_2B,...
        'a_symm_1B', a_symm_1B, 'moment_type', moment_type,...
        'resolution_error', resolution_error, 'force_pure_phases', force_pure_phases,...
        'spectral_algorithm_version', spectral_algorithm_version);
    disp("n = " + n);
    toc
    run_time = toc;
end
%% Plot the errors
energy_of_a = sum(abs(a_symm_1B) .^ 2);
SNR = energy_of_a ./ ((2*B + 1) * Q * sigma.^2);
error_squared_spectral_relative = error_squared_spectral ./ energy_of_a;
error_squared_FM_relative = error_squared_FM ./ energy_of_a;

prctile_region = 20; % 50 +- prctile_region percentile
asymp_start_percent = 80;

% Get asymptote value:
min_err_spectral = ...
    asymptoticSpectralError(Q, B, a_symm_1B, rho_hat_symm_2B, resolution_error);
min_err_spectral_relative = min_err_spectral ./ energy_of_a;

fig = plotMRAError('variable', N_vec, 'num_rep', num_rep_N, 'mean_or_median', 'median',...
'error_spectral', error_squared_spectral_relative, 'error_FM', error_squared_FM_relative,...
'labels', {"n", "Relative Squared Error"}, 'title', '', ...
'std_factor_or_prctile_region', prctile_region, 'is_x_math', true, ...
'y_asymptote_start_percent', asymp_start_percent, 'y_asymptote_val', min_err_spectral_relative);


%% Save the plot:

% Push test_number up by 1, make directory of test
dir_path = makeNewTestDir();
notes = "sigma = 0.1, made sure debias is correct";
figure_name = "changing_N_uniform_PS";
saveFigureToAllFormats(fig, figure_name, dir_path, notes);
sizeThresholdKB = 2000; %CHECK BEFORE LARGE RUN!!!!
saveNumericVariablesBelowThreshold(dir_path, "all_numeric_variables", sizeThresholdKB);
