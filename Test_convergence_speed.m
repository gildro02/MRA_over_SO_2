close all
%% Define Variables
N = 1e5;
B = 10;
Q = 2;
size_image = 51;
sigma = 20;
f = 0.1;
image_folder='Flower_Images';
image_name='Oxalis_tetraphylla_flower.jpg';
UPS = 0;
[a_symm_1B,image,Phi_ns_mat] = ...
    Generate_Picture_cut(['.\' image_folder '\' image_name],B,Q,size_image,UPS);

%% Relevant frequencies
[freq_file_name] = saveFrequencies(Q,B);
load(freq_file_name);

%% Rho Definition

temp = [0.00195583374369666 + 0.00530869159003379i;0.00195583374369666 + 0.00474988194897760i;0.00195583374369666 + 0.00419107230792141i;0.00195583374369666 + 0.00363226266686522i;0.00195583374369666 + 0.00307345302580903i;0.00195583374369666 + 0.00251464338475285i;0.00195583374369666 + 0.00195583374369666i;0.00195583374369666 + 0.00139702410264047i;0.00195583374369666 + 0.000838214461584282i;0.00195583374369666 + 0.000279404820528094i;0.00195583374369666 - 0.000279404820528094i;0.00195583374369666 - 0.000838214461584282i;0.00195583374369666 - 0.00139702410264047i;0.00195583374369666 - 0.00195583374369666i;0.00195583374369666 - 0.00251464338475285i;0.00195583374369666 - 0.00307345302580903i;0.00195583374369666 - 0.00363226266686522i;0.00195583374369666 - 0.00419107230792141i;0.00195583374369666 - 0.00474988194897760i;0.00195583374369666 - 0.00530869159003379i]...
    .*exp(1i*(f*rand(2*B,1)).^0.2)*3;

%temp = 5* (1: 2*B).';
t = -pi:1e-4:pi;
[rho_hat_symm_2B, rho_func] = GenerateDistribution(temp,t,Q,B);

rho_hat_symm_1B = rho_hat_symm_2B(-B<=k_symm_2B & k_symm_2B<=B);
rho_hat_half_2B = rho_hat_symm_2B(0<=k_symm_2B);
rho_hat_half_1B = rho_hat_symm_2B(0<=k_symm_2B & k_symm_2B<=B);
varrho = repelem(rho_hat_symm_1B,Q);

%% Generate Data

rotations=randsample(t,N,true,abs(rho_func(t))).';

a_symm_1B_mat = mat(a_symm_1B,Q,B);
a_half_1B_mat = a_symm_1B_mat(:,k_symm_1B>=0);
a_half_1B_vec = vec(a_half_1B_mat);

a_half_1B_mat_rot = a_half_1B_mat.*exp(-1i*k_half_1B_Q_mat.*permute(rotations,[3,2,1]))...
    +sigma*randn([size(a_half_1B_mat) N]);
a_half_1B_mat_rot_conj = conj(flip(reshape(...
    a_half_1B_mat_rot(repmat(k_half_1B_Q_mat~=0,[1,1,N]))...
    ,[Q,B,N]),2));
a_symm_1B_mat_rot = cat(2,a_half_1B_mat_rot_conj,a_half_1B_mat_rot);
a_symm_1B_rot = reshape(a_symm_1B_mat_rot,[Q*(2*B+1),N]);

%% True Moments
%Bias occurs only for q1=q2,k1=+-k2
Logical_Bias_Mat = zeros(length(a_symm_1B));
Logical_Bias_Mat(k_symm_1B_Q_vec-k_symm_1B_Q_vec.' == 0&...
    q_full_1B_vec-q_full_1B_vec.' == 0) = 1; %q1=q2,k1=k2
Logical_Bias_Mat(k_symm_1B_Q_vec+k_symm_1B_Q_vec.'==0&...
    q_full_1B_vec-q_full_1B_vec.' == 0) = 1; %q1=q2,k1=-k2

M_1_true = a_symm_1B.*varrho.*(2*pi);

T = BTTB(repmat(rho_hat_half_2B.', Q, 1));
M_2_true = (2*pi)*diag(a_symm_1B)*T*diag(conj(a_symm_1B))...
    +(sigma^2)*Logical_Bias_Mat; %'=complex conjugate transpose

%% Calculate Errors
M_1_Error = zeros(N,1);
M_2_Error = zeros(N,1);

M_1_Running = zeros(length(a_symm_1B), 1);
M_2_Running = zeros(length(a_symm_1B), length(a_symm_1B));
for n = 1:N
    M_1_Running = ((n-1) * M_1_Running + a_symm_1B_rot(:,n)) / n;
    M_2_Running = ((n-1) * M_2_Running + a_symm_1B_rot(:,n) * a_symm_1B_rot(:,n)') / n;
    %x*x' (not x*x.')!!!!
    
    M_1_Error(n) = norm(M_1_Running - M_1_true, "fro");
    M_2_Error(n) = norm(M_2_Running - M_2_true, "fro");
    
    if mod(n, 500) == 0
        done_percent = round(100 * n/N, 2);
        disp(done_percent + "% DONE");
    end
end

M_1_Error_Relative = M_1_Error ./ norm(a_symm_1B, "fro");
M_2_Error_Relative = M_2_Error ./ norm(a_symm_1B * a_symm_1B', "fro");

%% Plots
SNR = norm(a_symm_1B).^2 ./ (sigma.^2*(2*B+1)*Q);
figure
loglog(1:N, M_1_Error_Relative...
    , "LineWidth", 0.5);
hold on
loglog(1:N, M_2_Error_Relative...
    , "LineWidth", 0.5);
loglog(1:N,(1:N).^(-1/2)...
    , "LineWidth", 1.5); %slope -1/2 line

legend("Relative Error of M_1", "Relative Error of M_2", "Slope -1/2"...
    , "Fontsize", 15)
xlabel("N = Number of Realizations"...
    , "Fontsize", 15)
ylabel("Relative Error"...
    , "Fontsize", 15)
title("Asymptotic Behavior of Moment Error, SNR=" + SNR...
    , "Fontsize", 20)