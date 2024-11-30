%% Define variables
N = 1e5;
B = 10;
Q = 2;
t = -pi:1e-4:pi;
resolution_error = 1e-4;
moment_type = 'empirical';
%moment_type = 'analytical';
force_pure_phases = false;

num_unique_sigma = 90;
num_rep_sigma = 80;
sigma_vec_reduced = logspace(-3, 1, num_unique_sigma).';
sigma_vec = repelem(sigma_vec_reduced, num_rep_sigma);

%% Generate distribution:
f = 0;
temp = [0.00195583374369666 + 0.00530869159003379i;0.00195583374369666 + 0.00474988194897760i;0.00195583374369666 + 0.00419107230792141i;0.00195583374369666 + 0.00363226266686522i;0.00195583374369666 + 0.00307345302580903i;0.00195583374369666 + 0.00251464338475285i;0.00195583374369666 + 0.00195583374369666i;0.00195583374369666 + 0.00139702410264047i;0.00195583374369666 + 0.000838214461584282i;0.00195583374369666 + 0.000279404820528094i;0.00195583374369666 - 0.000279404820528094i;0.00195583374369666 - 0.000838214461584282i;0.00195583374369666 - 0.00139702410264047i;0.00195583374369666 - 0.00195583374369666i;0.00195583374369666 - 0.00251464338475285i;0.00195583374369666 - 0.00307345302580903i;0.00195583374369666 - 0.00363226266686522i;0.00195583374369666 - 0.00419107230792141i;0.00195583374369666 - 0.00474988194897760i;0.00195583374369666 - 0.00530869159003379i]...
    .* exp(1i * (f * rand(2*B, 1)) .^ 0.2) * 3;
%temp = [0.00195583374369666 + 0.00530869159003379i;0.00195583374369666 + 0.00474988194897760i;0.00195583374369666 + 0.00419107230792141i;0.00195583374369666 + 0.00363226266686522i;0.00195583374369666 + 0.00307345302580903i;0.00195583374369666 + 0.00251464338475285i;0.00195583374369666 + 0.00195583374369666i;0.00195583374369666 + 0.00139702410264047i;0.00195583374369666 + 0.000838214461584282i;0.00195583374369666 + 0.000279404820528094i;0.00195583374369666 - 0.000279404820528094i;0.00195583374369666 - 0.000838214461584282i;0.00195583374369666 - 0.00139702410264047i;0.00195583374369666 - 0.00195583374369666i;0.00195583374369666 - 0.00251464338475285i;0.00195583374369666 - 0.00307345302580903i;0.00195583374369666 - 0.00363226266686522i;0.00195583374369666 - 0.00419107230792141i;0.00195583374369666 - 0.00474988194897760i;0.00195583374369666 - 0.00530869159003379i];
[rho_hat_symm_2B, rho_func] = GenerateDistribution(temp, t, Q, B);

%% Make image:
size_image = 51;
image_folder = 'Flower_Images';
image_name = 'Oxalis_tetraphylla_flower.jpg';
[a_symm_1B, image, Phi_ns_mat] = ...
    Generate_Picture_cut(['.\' image_folder '\' image_name], B, Q, size_image, 0);

%% Initialize error vectors:
error_squared_FM = zeros(length(sigma_vec), 1);
error_squared_spectral = zeros(length(sigma_vec), 1);
tic

%% Run the experiment
for n = 1:length(sigma_vec)
    [error_squared_FM(n), coeff_FM, error_squared_spectral(n), coeff_spectral] = ...
        runNumericalExperiment('B', B, 'Q', Q, 'N', N, 'sigma', sigma_vec(n), 'rho_func',...
        rho_func, 't', t, 'rho_hat_symm_2B', rho_hat_symm_2B,...
        'a_symm_1B', a_symm_1B, 'moment_type', moment_type,...
        'resolution_error', resolution_error, 'force_pure_phases', force_pure_phases);
    disp("n = " + n);
    toc
end

%% Plot the errors

% Take the mean for each specific sigma:
error_squared_spectral_mean = mean(reshape(error_squared_spectral, [num_rep_sigma, num_unique_sigma]), 1).';
error_squared_FM_mean = mean(reshape(error_squared_FM, [num_rep_sigma, num_unique_sigma]), 1).';

% Define the SNR:
energy_of_a = sum(abs(a_symm_1B) .^ 2);
SNR = energy_of_a ./ ((2*B + 1) * Q * sigma_vec_reduced .^ 2);
fig = figure;
% plot Relative error:
loglog(SNR, error_squared_spectral_mean ./ energy_of_a)
hold on
loglog(SNR, error_squared_FM_mean./ energy_of_a)
xlabel("$\textnormal{SNR} = \frac{E_a}{(2B+1) Q \sigma^2}$", "fontsize", 15, "interpreter", "latex")
ylabel("$\textnormal{Relative Squared Error}$", "fontsize", 15, "interpreter", "latex")
axis tight
legend("Spectral Algorithm", "Frequancy Marching algorithm");

%% Save the plot:
% Push test_number up by 1, make directory of test
dir_path = makeNewTestDir();
notes = "N = 1e6;";
figure_name = "func_of_sigma.fig";
saveFigureAndNotes(fig, figure_name, dir_path, notes);
sizeThresholdKB = 1000;
saveNumericVariablesBelowThreshold(dir_path, "all_numeric_variables", sizeThresholdKB);



