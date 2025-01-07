%% Define variables
B = 10;
Q = 2;
t = -pi:1e-4:pi;
resolution_error_for_bound_algorithm = 1; %DOESNT MATTER
resolution_error = 5e-6; %FOR ACTUAL ERROR CALCULATIONS
moment_type = 'analytical'; %MUST BE ANALYTICAL
force_pure_phases = false;
sigma = 0; %DOESENT MATTER
N = 1e3; %DOESNT MATTER

%% Make image:
size_image = 51;
image_folder = 'Flower_Images';
image_name = 'Oxalis_tetraphylla_flower.jpg';
[a_symm_1B, image, Phi_ns_mat] = ...
    Generate_Picture_cut(['.\' image_folder '\' image_name], B, Q, size_image, 0);

%% Generate distributions & Calculate bound and error:
f = 0; %circulant
num_f = 200;
f_vec = logspace(3, 1, num_f).';
%f_vec=logspace(log10(500*B), log10(10*B), num_f).';

dist_from_circ = zeros(length(f_vec), 1);
bound = zeros(length(f_vec), 1);
error_squared_FM = zeros(length(f_vec), 1);
error_squared_spectral = zeros(length(f_vec), 1);
error_squared_spectral_unrestricted = zeros(length(f_vec), 1);
M_1_error = zeros(length(f_vec), 1);
M_2_error = zeros(length(f_vec), 1);
amp = 10;

% Generate Circulant distribution to perturb:
[rho_hat_symm_2B_circ, rho_func_circ] =...
    generateCirculantDistribution(GenerateDistribution(rand(2*B, 1), t, Q, B), t, Q, B);
tic
for n=1:length(f_vec)
    
    % Generate distribution:
    %temp = [0.00195583374369666 + 0.00530869159003379i;0.00195583374369666 + 0.00474988194897760i;0.00195583374369666 + 0.00419107230792141i;0.00195583374369666 + 0.00363226266686522i;0.00195583374369666 + 0.00307345302580903i;0.00195583374369666 + 0.00251464338475285i;0.00195583374369666 + 0.00195583374369666i;0.00195583374369666 + 0.00139702410264047i;0.00195583374369666 + 0.000838214461584282i;0.00195583374369666 + 0.000279404820528094i;0.00195583374369666 - 0.000279404820528094i;0.00195583374369666 - 0.000838214461584282i;0.00195583374369666 - 0.00139702410264047i;0.00195583374369666 - 0.00195583374369666i;0.00195583374369666 - 0.00251464338475285i;0.00195583374369666 - 0.00307345302580903i;0.00195583374369666 - 0.00363226266686522i;0.00195583374369666 - 0.00419107230792141i;0.00195583374369666 - 0.00474988194897760i;0.00195583374369666 - 0.00530869159003379i]...
    %    .* exp(1i * rand(2*B, 1) .^ 0.2 ./ f_vec(n)) * amp;
    
    temp = rho_hat_symm_2B_circ(2*B + 2:end) .* exp(1i * rand(2*B, 1) .^ 0.2 ./ f_vec(n)) * amp;
    [rho_hat_symm_2B, rho_func] = GenerateDistribution(temp, t, Q, B);
    dist_from_circ(n) = distanceFromCirculant(Q, B, rho_hat_symm_2B);
    
    % Calculate bound:
    bound(n) = boundForSpectralAlgorithm_NewVersion('B', B, 'Q', Q, 'N', N, 'sigma', sigma, 'rho_func',...
        rho_func, 't', t, 'rho_hat_symm_2B', rho_hat_symm_2B,...
        'a_symm_1B', a_symm_1B, 'moment_type', moment_type,...
        'resolution_error', resolution_error_for_bound_algorithm, 'force_pure_phases', force_pure_phases);
    
    % Calculate errors for analytical Moments:
    [~, ~, error_squared_spectral(n), ~, M_1_error(n), M_2_error(n)] = ...
        runNumericalExperiment('B', B, 'Q', Q, 'N', N, 'sigma', sigma, 'rho_func',...
        rho_func, 't', t, 'rho_hat_symm_2B', rho_hat_symm_2B,...
        'a_symm_1B', a_symm_1B, 'moment_type', moment_type,...
        'resolution_error', resolution_error, 'force_pure_phases', force_pure_phases,...
        'spectral_error_func', 'restricted', 'only_one_algorithm', 'Spectral');
    
    % Calculate errors with unrestricted error function:
    [~, ~, error_squared_spectral_unrestricted(n), ~, M_1_error(n), M_2_error(n)] = ...
        runNumericalExperiment('B', B, 'Q', Q, 'N', N, 'sigma', sigma, 'rho_func',...
        rho_func, 't', t, 'rho_hat_symm_2B', rho_hat_symm_2B,...
        'a_symm_1B', a_symm_1B, 'moment_type', moment_type,...
        'resolution_error', resolution_error, 'force_pure_phases', force_pure_phases,...
        'spectral_error_func', 'unrestricted', 'only_one_algorithm', 'Spectral');

    disp("n = " + n);
    toc
end

%% Plot
energy_of_a = sum(abs(a_symm_1B).^2);
fig = figure;
loglog(dist_from_circ, bound./energy_of_a, "-", "LineWidth", 1.5);
hold on
loglog(dist_from_circ, error_squared_spectral./energy_of_a, ".", "MarkerSize", 10);
loglog(dist_from_circ, error_squared_spectral_unrestricted./energy_of_a, ".", "MarkerSize", 10);
legend("Bound", "Steps of $\frac{1}{2B+1}$", "Steps of $10^{-6}$",...
    "FontSize", 15, "interpreter", "latex");
xlabel("$||T-C||^2$", "FontSize", 15, "Interpreter", "latex")
ylabel("Relative Squared Error", "Interpreter", "latex", "FontSize", 15)

set(gca, 'xdir', 'reverse')
axis tight
grid on


%% Save the plot:
%{
% Push test_number up by 1, make directory of test
dir_path = makeNewTestDir();
figure_name = 'bound_test_with_both_errors';
saveFigureAndNotes(fig, [figure_name, '.fig'], dir_path);
saveFigureToEPS(fig, [figure_name, '.eps'], dir_path);
sizeThresholdKB = 1000;
saveNumericVariablesBelowThreshold(dir_path, "all_numeric_variables", sizeThresholdKB);
%}