B=10;
Q=2;
%% Relevant frequencies
k_symm_1B=-B:B;
k_half_1B=0:B;
k_symm_2B=-2*B:2*B;
k_half_2B=0:2*B;

k_symm_1B_Q_mat=repmat(k_symm_1B,[Q,1]);
k_half_1B_Q_mat=repmat(k_half_1B,[Q,1]);
k_symm_2B_Q_mat=repmat(k_symm_2B,[Q,1]);
k_half_2B_Q_mat=repmat(k_half_2B,[Q,1]);

k_symm_1B_Q_vec=vec(k_symm_1B_Q_mat);
k_half_1B_Q_vec=vec(k_half_1B_Q_mat);
k_symm_2B_Q_vec=vec(k_symm_2B_Q_mat);
k_half_2B_Q_vec=vec(k_half_2B_Q_mat);

q_full_1B_mat=repmat((0:Q-1).',[1,2*B+1]);
q_full_1B_vec=vec(q_full_1B_mat);
q_1B=0:Q-1;
%% Definition of Picture
%{
isUniformPowerSpectrum=[0,1];
images=cell(length(isUniformPowerSpectrum),1);
a_symm_1B_cell=cell(length(isUniformPowerSpectrum),1);
%}
%[a_symm_1B,image]=Generate_Picture_cut('flower-1.jpg',B,Q,false);
[a_symm_1B,image]=Generate_Picture_cut('flower-1.jpg',B,Q,true);
a_symm_1B_mat=mat(a_symm_1B,Q,B);
a_half_1B_mat=a_symm_1B_mat(:,k_symm_1B>=0);
a_half_1B_vec=vec(a_half_1B_mat);
%% Rho Definition
%temp=[0.00195583374369666 + 0.00530869159003379i;0.00195583374369666 + 0.00474988194897760i;0.00195583374369666 + 0.00419107230792141i;0.00195583374369666 + 0.00363226266686522i;0.00195583374369666 + 0.00307345302580903i;0.00195583374369666 + 0.00251464338475285i;0.00195583374369666 + 0.00195583374369666i;0.00195583374369666 + 0.00139702410264047i;0.00195583374369666 + 0.000838214461584282i;0.00195583374369666 + 0.000279404820528094i;0.00195583374369666 - 0.000279404820528094i;0.00195583374369666 - 0.000838214461584282i;0.00195583374369666 - 0.00139702410264047i;0.00195583374369666 - 0.00195583374369666i;0.00195583374369666 - 0.00251464338475285i;0.00195583374369666 - 0.00307345302580903i;0.00195583374369666 - 0.00363226266686522i;0.00195583374369666 - 0.00419107230792141i;0.00195583374369666 - 0.00474988194897760i;0.00195583374369666 - 0.00530869159003379i]...
%    .*exp(1i*(f_vec(f_index)*rand(2*B,1)).^0.2)*3;
temp=[0.00195583374369666 + 0.00530869159003379i;0.00195583374369666 + 0.00474988194897760i;0.00195583374369666 + 0.00419107230792141i;0.00195583374369666 + 0.00363226266686522i;0.00195583374369666 + 0.00307345302580903i;0.00195583374369666 + 0.00251464338475285i;0.00195583374369666 + 0.00195583374369666i;0.00195583374369666 + 0.00139702410264047i;0.00195583374369666 + 0.000838214461584282i;0.00195583374369666 + 0.000279404820528094i;0.00195583374369666 - 0.000279404820528094i;0.00195583374369666 - 0.000838214461584282i;0.00195583374369666 - 0.00139702410264047i;0.00195583374369666 - 0.00195583374369666i;0.00195583374369666 - 0.00251464338475285i;0.00195583374369666 - 0.00307345302580903i;0.00195583374369666 - 0.00363226266686522i;0.00195583374369666 - 0.00419107230792141i;0.00195583374369666 - 0.00474988194897760i;0.00195583374369666 - 0.00530869159003379i]...
    .*exp(1i*(rand(2*B,1)).^0.2.*0)*3;
%temp=[0.00195583374369666 + 0.00530869159003379i;0.00195583374369666 + 0.00474988194897760i;0.00195583374369666 + 0.00419107230792141i;0.00195583374369666 + 0.00363226266686522i;0.00195583374369666 + 0.00307345302580903i;0.00195583374369666 + 0.00251464338475285i;0.00195583374369666 + 0.00195583374369666i;0.00195583374369666 + 0.00139702410264047i;0.00195583374369666 + 0.000838214461584282i;0.00195583374369666 + 0.000279404820528094i;0.00195583374369666 - 0.000279404820528094i;0.00195583374369666 - 0.000838214461584282i;0.00195583374369666 - 0.00139702410264047i;0.00195583374369666 - 0.00195583374369666i;0.00195583374369666 - 0.00251464338475285i;0.00195583374369666 - 0.00307345302580903i;0.00195583374369666 - 0.00363226266686522i;0.00195583374369666 - 0.00419107230792141i;0.00195583374369666 - 0.00474988194897760i;0.00195583374369666 - 0.00530869159003379i]...
%.*exp(1i*(0.8*ones(2*B,1)).^0.3./f_vec(f_index).^2)*3;
%temp=[0.00195583374369666 + 0.00530869159003379i;0.00195583374369666 + 0.00474988194897760i;0.00195583374369666 + 0.00419107230792141i;0.00195583374369666 + 0.00363226266686522i;0.00195583374369666 + 0.00307345302580903i;0.00195583374369666 + 0.00251464338475285i;0.00195583374369666 + 0.00195583374369666i;0.00195583374369666 + 0.00139702410264047i;0.00195583374369666 + 0.000838214461584282i;0.00195583374369666 + 0.000279404820528094i;0.00195583374369666 - 0.000279404820528094i;0.00195583374369666 - 0.000838214461584282i;0.00195583374369666 - 0.00139702410264047i;0.00195583374369666 - 0.00195583374369666i;0.00195583374369666 - 0.00251464338475285i;0.00195583374369666 - 0.00307345302580903i;0.00195583374369666 - 0.00363226266686522i;0.00195583374369666 - 0.00419107230792141i;0.00195583374369666 - 0.00474988194897760i;0.00195583374369666 - 0.00530869159003379i]...
%    +f(index)/(5e4*B)*(rand(2*B,1)+1i*rand(2*B,1));
rho_hat_symm_2B_non_normalized =[conj(flip(temp)); 1/(2*pi); temp];
rho_hat_symm_2B_temp=(1/(2*pi)) * (1/abs(rho_hat_symm_2B_non_normalized(k_symm_2B==0)))...
    * rho_hat_symm_2B_non_normalized;
rho_func_temp = @(x) real(rho_hat_symm_2B_temp.'*exp(1i*x'*k_symm_2B).');

t=-pi:1e-4:pi;
min_rho_func=min(rho_func_temp(t));
rho_hat_symm_2B =rho_hat_symm_2B_temp;
rho_hat_symm_2B(k_symm_2B==0)=rho_hat_symm_2B(k_symm_2B==0)-min(min_rho_func,0); %make >=0
rho_hat_symm_2B=rho_hat_symm_2B.*(1/(2*pi))*(1/rho_hat_symm_2B(k_symm_2B==0)); %renormalize
rho_func = @(x) real(rho_hat_symm_2B.'*exp(1i*x.'*k_symm_2B).'); %the final rho_func

rho_hat_symm_1B=rho_hat_symm_2B(-B<=k_symm_2B & k_symm_2B<=B);
rho_hat_half_2B=rho_hat_symm_2B(0<=k_symm_2B);
rho_hat_half_1B=rho_hat_symm_2B(0<=k_symm_2B & k_symm_2B<=B);

varrho=repelem(rho_hat_symm_1B,Q);
%% Generate Data
N=1e5; %num of rotated photos;

dt=t(2)-t(1); %Quantization of [0,2pi].
rotations=randsample(t,N,true,abs(rho_func(t))).';
%absolute value to avoid numeric negatives, rho should always be positive.
sigma=1; %noise level


a_half_1B_mat_rot=a_half_1B_mat.*exp(-1i*k_half_1B_Q_mat.*permute(rotations,[3,2,1]))...
    +sigma*randn(size(a_half_1B_mat));
a_half_1B_mat_rot_conj=conj(flip(reshape(...
    a_half_1B_mat_rot(repmat(k_half_1B_Q_mat~=0,[1,1,N]))...
    ,[Q,B,N]),2));
a_symm_1B_mat_rot=cat(2,a_half_1B_mat_rot_conj,a_half_1B_mat_rot);
a_symm_1B_rot=reshape(a_symm_1B_mat_rot,[Q*(2*B+1),N]);

%a_half_1B_rot=a_half_1B_vec.*exp(1i.*k_half_1B_Q_vec.*rotations.');
%a_symm_1B_rot=[conj(reshape(flip(reshape(a_half_1B_rot(k_half_1B_Q_vec~=0),[Q,B,N]),2),[Q*B,N])) ; a_half_1B_rot];
%% Calculate Moments (M1 and M2)

%first moment
M_1_true=a_symm_1B.*varrho.*(2*pi);
M_1_emp=vec(mean(a_symm_1B_mat_rot,3));

%second moment
T=BTTB(repmat(rho_hat_half_2B.',Q,1));

%Bias occurs only for q1=q2,k1=+-k2
Logical_Bias_Mat=zeros(length(a_symm_1B));
Logical_Bias_Mat(k_symm_1B_Q_vec-k_symm_1B_Q_vec.'==0&...
    q_full_1B_vec-q_full_1B_vec.'==0) = 1; %q1=q2,k1=k2
Logical_Bias_Mat(k_symm_1B_Q_vec+k_symm_1B_Q_vec.'==0&...
    q_full_1B_vec-q_full_1B_vec.'==0) = 1; %q1=q2,k1=-k2

M_2_true=(2*pi)*diag(a_symm_1B)*T*diag(conj(a_symm_1B))+(sigma^2)*Logical_Bias_Mat; %'=complex conjugate transpose
M_2_emp=(1/N)*(a_symm_1B_rot*a_symm_1B_rot');

M_2_true_wave=M_2_true-(sigma^2)*Logical_Bias_Mat;
M_2_emp_wave=M_2_emp-(sigma^2)*Logical_Bias_Mat;

%moment_type="emp";
moment_type="true";
if strcmp(moment_type,"emp")==1
    M_1=M_1_emp;
    M_2=M_2_emp;
    M_2_wave=M_2_emp_wave;
end
if strcmp(moment_type,"true")==1
    M_1=M_1_true;
    M_2=M_2_true;
    M_2_wave=M_2_true_wave;
end
%% Frequency Marching Algorithm
S_mathcal=2*pi*diag(1./M_1)*M_2_wave*diag(1./M_1)';
weights_q=cell(2*B+1,2*B+1);
for k1=-B:B
    for k2=-B:B
        weight=exp(-q_1B).'*exp(-q_1B); %low q's are worth most
        weights_q{k_symm_1B==k1,k_symm_1B==k2}=weight./sum(weight,"all"); %normalize
    end
end
S_mathcal_tilde=zeros(2*B+1,2*B+1); %condensed S
for k1=-B:B
    for k2=-B:B
        S_mathcal_tilde(k_symm_1B==k1,k_symm_1B==k2)...
            =sum(weights_q{k_symm_1B==k1,k_symm_1B==k2}...
            .*S_mathcal(k_symm_1B_Q_vec==k1,k_symm_1B_Q_vec==k2),"all");
    end
end
weights_k=cell(1,B); %max number of approximations x max freq-1
for k=2:B
    weights_k{k}=exp(-abs(ceil(k/2)-(1:k-1).')); %lorentzian of length n centered at ceil(n/2)
    weights_k{k}=weights_k{k}./sum(weights_k{k}); %normalize to sum 1.
end
%%
rho_approx_half_1B=zeros(length(k_half_1B),1);
rho_approx_tilde_half_1B= zeros(length(k_half_1B),1);

rho_abs_half_1B = sqrt(abs(1./(2*pi*diag(S_mathcal_tilde(k_symm_1B>=0,k_symm_1B>=0)))));
rho_approx_half_1B(k_half_1B==0)=1/(2*pi);
%arbitrarily pick phase of 1'st coeff to be 1:
rho_approx_half_1B(k_half_1B==1) = rho_abs_half_1B(k_half_1B==1);

for k=2:B
    rho_approx_tilde_half_1B(k_half_1B==k) = sum(weights_k{k}.*flip(rho_approx_half_1B(1<=k_half_1B & k_half_1B <= k-1))./...
        conj(rho_approx_half_1B(1<=k_half_1B & k_half_1B <= k-1))./...
        S_mathcal_tilde(k_symm_1B==k, 1<=k_symm_1B & k_symm_1B <= k-1).');
    
    rho_approx_half_1B(k_half_1B==k)=rho_abs_half_1B(k_half_1B==k).*...
        rho_approx_tilde_half_1B(k_half_1B==k)./abs(rho_approx_tilde_half_1B(k_half_1B==k));
end
rho_approx_symm_1B= [conj(flip(rho_approx_half_1B(2:end))); rho_approx_half_1B];
varrho_approx=repelem(rho_approx_symm_1B,Q);
%[err_squared_rho,theta_min]=circ_error_continuous_unrestricted(rho_approx_symm_1B,rho_hat_symm_1B,B)

a_approx_symm_1B=M_1./(2*pi.*varrho_approx);
[err_squared_a_FM,theta_min]=circ_error_continuous_unrestricted_2D(a_approx_symm_1B,a_symm_1B,Q,B);

