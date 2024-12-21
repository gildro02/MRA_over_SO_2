%% Define variables
B = 10;
Q = 2;
t = -pi:1e-4:pi;
resolution_error = 1e-4;
moment_type = 'empirical';
%moment_type = 'analytical';
force_pure_phases = false;
sigma = 0;
num_unique_N = 75;
num_rep_N = 1000;
N_vec_reduced = ceil(logspace(0, 6, num_unique_N).');
N_vec = repelem(N_vec_reduced, num_rep_N);

%% Generate distribution:
f = 0; %circulant
amp = 1;
temp = [0.00195583374369666 + 0.00530869159003379i;0.00195583374369666 + 0.00474988194897760i;0.00195583374369666 + 0.00419107230792141i;0.00195583374369666 + 0.00363226266686522i;0.00195583374369666 + 0.00307345302580903i;0.00195583374369666 + 0.00251464338475285i;0.00195583374369666 + 0.00195583374369666i;0.00195583374369666 + 0.00139702410264047i;0.00195583374369666 + 0.000838214461584282i;0.00195583374369666 + 0.000279404820528094i;0.00195583374369666 - 0.000279404820528094i;0.00195583374369666 - 0.000838214461584282i;0.00195583374369666 - 0.00139702410264047i;0.00195583374369666 - 0.00195583374369666i;0.00195583374369666 - 0.00251464338475285i;0.00195583374369666 - 0.00307345302580903i;0.00195583374369666 - 0.00363226266686522i;0.00195583374369666 - 0.00419107230792141i;0.00195583374369666 - 0.00474988194897760i;0.00195583374369666 - 0.00530869159003379i]...
    .* exp(1i * (f * rand(2*B, 1)) .^ 0.2) * amp;
%temp = [0.00195583374369666 + 0.00530869159003379i;0.00195583374369666 + 0.00474988194897760i;0.00195583374369666 + 0.00419107230792141i;0.00195583374369666 + 0.00363226266686522i;0.00195583374369666 + 0.00307345302580903i;0.00195583374369666 + 0.00251464338475285i;0.00195583374369666 + 0.00195583374369666i;0.00195583374369666 + 0.00139702410264047i;0.00195583374369666 + 0.000838214461584282i;0.00195583374369666 + 0.000279404820528094i;0.00195583374369666 - 0.000279404820528094i;0.00195583374369666 - 0.000838214461584282i;0.00195583374369666 - 0.00139702410264047i;0.00195583374369666 - 0.00195583374369666i;0.00195583374369666 - 0.00251464338475285i;0.00195583374369666 - 0.00307345302580903i;0.00195583374369666 - 0.00363226266686522i;0.00195583374369666 - 0.00419107230792141i;0.00195583374369666 - 0.00474988194897760i;0.00195583374369666 - 0.00530869159003379i];
[rho_hat_symm_2B, rho_func] = GenerateDistribution(temp, t, Q, B);
dist_from_circ = distanceFromCirculant(Q, B, rho_hat_symm_2B);

%% Make image:
size_image = 51;
image_folder = 'Flower_Images';
image_name = 'Oxalis_tetraphylla_flower.jpg';
%image_name = 'flower-circle-border.jpg';
%image_name = 'nonzero_coefficients_offdiagonal.jpg';
[a_symm_1B, image, Phi_ns_mat] = ...
    Generate_Picture_cut(['.\' image_folder '\' image_name], B, Q, size_image, 0);
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
        'resolution_error', resolution_error, 'force_pure_phases', force_pure_phases);
    disp("n = " + n);
    toc
end
%% Plot the errors
energy_of_a = sum(abs(a_symm_1B) .^ 2);
SNR = energy_of_a ./ ((2*B + 1) * Q * sigma.^2);
error_squared_spectral_relative = error_squared_spectral ./ energy_of_a;
error_squared_FM_relative = error_squared_FM ./ energy_of_a;
fig = plotErrorAsFunctionOfX(N_vec, "N", error_squared_spectral_relative,...
    error_squared_FM_relative, num_rep_N, dist_from_circ);

%% Save the plot:
% Push test_number up by 1, make directory of test
dir_path = makeNewTestDir();
notes = "sigma = 0";
figure_name = "infinite_snr.fig";
saveFigureAndNotes(fig, figure_name, dir_path, notes);
sizeThresholdKB = 2000; %CHECK BEFORE LARGE RUN!!!!
saveNumericVariablesBelowThreshold(dir_path, "all_numeric_variables", sizeThresholdKB);