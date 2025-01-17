%% Define variables
B = 10;
Q = 2;
t = -pi:1e-4:pi;
sigma = 50;
%% Generate distribution:
f = 1; %circulant
amp = 5;
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

%% Build struct of parameters:
params.B = B;
params.Q = Q;
params.t = t;
params.sigma = sigma;
params.rho_hat_symm_2B = rho_hat_symm_2B;
params.rho_func = rho_func;
params.a_symm_1B = a_symm_1B;

%% Get moments
num_N_unique = 10;
num_rep_N = 10;
N_vec_reduced = ceil(logspace(0,6,num_N_unique));
N_vec = repelem(N_vec_reduced, num_rep_N);
diff_vec_M1 = zeros(num_N_unique, 1);
diff_vec_M2 = zeros(num_N_unique, 1);

tic
for N_idx=1:length(N_vec)
    params.N = N_vec(N_idx);
    [M_1_emp, M_2_emp] = generateEmpiricalMoments(params);
    [M_1_analytic, M_2_analytic] = generateAnalyticalMoments(params);
    diff_vec_M1(N_idx) = norm(M_1_emp - M_1_analytic, "fro");
    diff_vec_M2(N_idx) = norm(M_2_emp - M_2_analytic, "fro");
    disp("n = " + N_idx);
    toc
end
%% Plot
%close all
figure
diff_vec_M1_mean = mean(reshape(diff_vec_M1,...
    [num_rep_N, num_N_unique]), 1);
diff_vec_M2_mean = mean(reshape(diff_vec_M2,...
    [num_rep_N, num_N_unique]), 1);
loglog(N_vec_reduced, diff_vec_M1_mean)
hold on
loglog(N_vec_reduced, diff_vec_M2_mean)
loglog(N_vec_reduced, N_vec_reduced .^ (-1/2))
title("sigma = " + sigma)
legend("M_1", "M_2", "line -1/2")