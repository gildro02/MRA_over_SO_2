B=10;
Q=2;

[freq_file_name]=saveFrequencies(Q,B);
load(freq_file_name);
UPS=1;
[a_symm_1B,image]=...
    Generate_Picture_cut([pwd '\Flower_Images\Oxalis_tetraphylla_flower.jpg'],B,Q,UPS);

a_symm_1B_mat=mat(a_symm_1B,Q,B);
a_half_1B_mat=a_symm_1B_mat(:,k_symm_1B>=0);
a_half_1B_vec=vec(a_half_1B_mat);

f_vec=0;
sigma=1e-1;
repetitions=5;
N_max=1e6;

%temp=[0.00195583374369666 + 0.00530869159003379i;0.00195583374369666 + 0.00474988194897760i;0.00195583374369666 + 0.00419107230792141i;0.00195583374369666 + 0.00363226266686522i;0.00195583374369666 + 0.00307345302580903i;0.00195583374369666 + 0.00251464338475285i;0.00195583374369666 + 0.00195583374369666i;0.00195583374369666 + 0.00139702410264047i;0.00195583374369666 + 0.000838214461584282i;0.00195583374369666 + 0.000279404820528094i;0.00195583374369666 - 0.000279404820528094i;0.00195583374369666 - 0.000838214461584282i;0.00195583374369666 - 0.00139702410264047i;0.00195583374369666 - 0.00195583374369666i;0.00195583374369666 - 0.00251464338475285i;0.00195583374369666 - 0.00307345302580903i;0.00195583374369666 - 0.00363226266686522i;0.00195583374369666 - 0.00419107230792141i;0.00195583374369666 - 0.00474988194897760i;0.00195583374369666 - 0.00530869159003379i]...
%    .*exp(1i*(f_vec*rand(2*B,1)).^0.2)*3;
temp = (1/(2*pi))*ones(2*B,1);
t=-pi:1e-4:pi;
[rho_hat_symm_2B, rho_func]=GenerateDistribution(temp,t,Q,B);


rho_hat_symm_1B=rho_hat_symm_2B(-B<=k_symm_2B & k_symm_2B<=B);
rho_hat_half_2B=rho_hat_symm_2B(0<=k_symm_2B);
rho_hat_half_1B=rho_hat_symm_2B(0<=k_symm_2B & k_symm_2B<=B);


T=BTTB(repmat(rho_hat_half_2B.',Q,1));
%Bias occurs only for q1=q2,k1=+-k2

Logical_Bias_Mat=zeros(length(a_symm_1B));
Logical_Bias_Mat(k_symm_1B_Q_vec-k_symm_1B_Q_vec.'==0&...
    q_full_1B_vec-q_full_1B_vec.'==0) = 1; %q1=q2,k1=k2
Logical_Bias_Mat(k_symm_1B_Q_vec+k_symm_1B_Q_vec.'==0&...
    q_full_1B_vec-q_full_1B_vec.'==0) = 1; %q1=q2,k1=-k2

%Logical_bias_Mat=eye(length(a_symm_1B))
M_2_true=(2*pi)*diag(a_symm_1B)*T*diag(conj(a_symm_1B))+(sigma^2)*Logical_Bias_Mat; %'=complex conjugate transpose
%M_2_true=(2*pi)*diag(a_symm_1B)*T*diag(conj(a_symm_1B))+(1/(2*pi))*(sigma^2)*Logical_Bias_Mat; %'=complex conjugate transpose

error_sq_vec=zeros(N_max,repetitions);
M_2_emp=0;

for rep=1:repetitions
    dt=t(2)-t(1); %Quantization of [0,2pi].
    rotations=randsample(t,N_max,true,abs(rho_func(t))).';
    %rotataions=zeros(N_max,1);
    %absolute value to avoid numeric negatives, rho should always be positive.


    a_half_1B_mat_rot=a_half_1B_mat.*exp(-1i*k_half_1B_Q_mat.*permute(rotations,[3,2,1]))...
        + sigma*randn([size(a_half_1B_mat) N_max]);
    %a_half_1B_mat_rot=a_half_1B_mat.*exp(-1i*k_half_1B_Q_mat.*permute(rotations,[3,2,1]))...
    %    +(sigma/sqrt(2*pi))*randn(size(a_half_1B_mat));

    a_half_1B_mat_rot_conj=conj(flip(reshape(...
        a_half_1B_mat_rot(repmat(k_half_1B_Q_mat~=0,[1,1,N_max]))...
        ,[Q,B,N_max]),2));
    a_symm_1B_mat_rot=cat(2,a_half_1B_mat_rot_conj,a_half_1B_mat_rot);
    a_symm_1B_rot=reshape(a_symm_1B_mat_rot,[Q*(2*B+1),N_max]);

    for N=1:N_max
        M_2_emp=(1/N)*((N-1)*M_2_emp+a_symm_1B_rot(:,N)*a_symm_1B_rot(:,N)');
        error_sq_vec(N,rep)=norm(M_2_emp-M_2_true,"fro").^2;
        if mod(N,1e4)==0
            disp("running at "+N+", on repetition "+rep);
        end
    end
end
error_sq_vec_avg=mean(error_sq_vec,2);
%disp(error_sq_vec);
figure
loglog(1:N_max,error_sq_vec_avg)

M_2_diff=M_2_emp-M_2_true;