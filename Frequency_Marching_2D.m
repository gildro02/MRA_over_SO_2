close all
addpath(genpath('basis'))
addpath(genpath('common'))
addpath(genpath('examples'))
addpath(genpath('fourier'))
addpath(genpath('projections'))
%% Initialyzing - comment after first run
tic
if ~exist('picture','var')
picture=imread('flower-1.jpg');
Full_Picture=im2double(rgb2gray(imresize(picture,51/1000)));
L=size(Full_Picture,1);
N=floor(L/2);
[x, y]=meshgrid(-N:N, -N:N);
r=sqrt(x.^2+y.^2);
Compressed_Picture=Full_Picture(r<=N);
Full_Picture(r>N)=0;

% Generate Basis:

[ Phi_ns, ang_freqs, rad_freqs, R_ns ]=Bessel_ns_v5( N );
Phi_ns_mat=[Phi_ns, conj(Phi_ns(:, ang_freqs~=0))];
coeff = Phi_ns_mat\Compressed_Picture;
full_ang_freqs=[ang_freqs;-ang_freqs(ang_freqs~=0)];
full_rad_freqs=[rad_freqs;rad_freqs(ang_freqs~=0)];
Compressed_Proj_Picture=Phi_ns_mat*coeff;
Proj_Picture=zeros(L);
Proj_Picture(r<=N)=Compressed_Proj_Picture;
Proj_Picture=real(Proj_Picture); %%picture we will use

ang_freqs_max=max(ang_freqs);
B=ang_freqs_max;
Q_vec=zeros(B+1,1); %number of radial frequencies per angular frequency
                                %from 0 to B
for index1=0:ang_freqs_max
    Q_vec(index1+1)=sum(ang_freqs==index1);
end
end
toc
%% Generate Data
B=ang_freqs_max;%20 %Band Limit
seed_vec=1./(1:B/2).^2;
[rho_coeff_symm,rho_func]=Generate_Rho_Squared(seed_vec);
rho_coeff_symm_indecies=-B:B;
rho_coeff=rho_coeff_symm(rho_coeff_symm_indecies>=0);
P=10000; %num of rotated photos;

dx=1e-5; %Quantization of [0,2pi].
rot_angles=randsample(0:dx:2*pi,P,true,dx*rho_func(0:dx:2*pi))'; %randomize from distribution[0,2pi].

sigma=1; %noise level

%calculate noise:
coeff_Noise_positive=sigma*randn(length(ang_freqs),P);
coeff_Noise=[coeff_Noise_positive; coeff_Noise_positive(ang_freqs~=0,:)];
%Noisy coefficients: original coeff * phases + real noise:
coeff_data=repmat(coeff,1,P).*exp(-1i*full_ang_freqs*(rot_angles')) + coeff_Noise;

%% Calculate Moments (M1 and M2)
%rho_coeff(k), but copied over the q dimension.
altered_rho_coeff=zeros(length(coeff),1);
for index1=-ang_freqs_max:ang_freqs_max
    temp=rho_coeff_symm(rho_coeff_symm_indecies==index1);
    if ~isempty(temp)
    altered_rho_coeff(full_ang_freqs==index1)=temp;
    end
end
%first moment
M1_anal=coeff.*altered_rho_coeff.*(2*pi); %*(2pi)
M1_emp=mean(coeff_data,2);

%second moment
ang_freq_index_Mat=full_ang_freqs-full_ang_freqs.'; %colomn minus row
rad_freq_index_Mat=full_rad_freqs-full_rad_freqs.';
rho_coeff_Mat=zeros(length(coeff));
%place correct values to rho_coeff_Mat
for index1=-ang_freqs_max:ang_freqs_max
    temp=rho_coeff_symm(rho_coeff_symm_indecies==index1);
    if ~isempty(temp)
    rho_coeff_Mat(ang_freq_index_Mat==index1)=temp;
    end
end
%Bias occurs only for q1=q2,k1=+-k2
Logical_Bias_Mat=zeros(length(coeff));
Logical_Bias_Mat((rad_freq_index_Mat==0)&ang_freq_index_Mat==0)=1;
Logical_Bias_Mat((rad_freq_index_Mat==0)&(full_ang_freqs+full_ang_freqs.'==0))=1;

M2_anal=coeff*coeff'.*rho_coeff_Mat.*(2*pi)+(sigma^2)*Logical_Bias_Mat; %'=complex conjugate transpose
M2_emp=(1/P)*(coeff_data*coeff_data');

%% Frequency Marching Algorithm

M2_wave_anal=2*pi*diag(1./M1_anal)*...
    (M2_anal-(sigma^2)*Logical_Bias_Mat)*...
    (diag(1./M1_anal))';
M2_wave_anal_positive_freq=M2_wave_anal(1:length(ang_freqs),1:length(ang_freqs));

M2_wave_emp=2*pi*diag(1./M1_emp)*...
    (M2_emp-(sigma^2)*Logical_Bias_Mat)*...
    (diag(1./M1_emp))';
M2_wave_emp_positive_freq=M2_wave_emp(1:length(ang_freqs),1:length(ang_freqs));

M2_wave_compact=zeros(ang_freqs_max+1,ang_freqs_max+1);
%weights=cell(
freq_compact=0:ang_freqs_max;
for index1=0:ang_freqs_max
    for index2=0:ang_freqs_max
        M2_wave_compact(freq_compact==index1,freq_compact==index2)=...
            1./mean(1./M2_wave_anal_positive_freq(ang_freqs==index1,ang_freqs==index2),"all");
    end
end

rho_coeff_approx=zeros(1,ang_freqs_max+1);
rho_coeff_approx(freq_compact==0)=1/(2*pi);
rho_coeff_approx(freq_compact==1)=1/sqrt(abs(...
    M2_wave_compact(freq_compact==1,freq_compact==1)/rho_coeff_approx(freq_compact==0)));
for n=2:ang_freqs_max
    rho_coeff_approx_raw=1./(M2_wave_compact(freq_compact==n,(1<=freq_compact)&(freq_compact<n))).*...
                 flip(rho_coeff_approx((1<=freq_compact)&(freq_compact<n)))./...
                 conj(rho_coeff_approx((1<=freq_compact)&(freq_compact<n)));
     mean_rho_coeff_approx_raw=mean(rho_coeff_approx_raw);
     rho_coeff_approx(freq_compact==n)=mean_rho_coeff_approx_raw;
end

rho_coeff_symm_approx=[flip(conj(rho_coeff_approx(2:end))), rho_coeff_approx];
altered_rho_coeff_approx=zeros(length(coeff),1);
for index1=-ang_freqs_max:ang_freqs_max
    temp=rho_coeff_symm_approx(rho_coeff_symm_indecies==index1);
    if ~isempty(temp)
    altered_rho_coeff_approx(full_ang_freqs==index1)=temp;
    end
end

coeff_approx=(M1_anal./altered_rho_coeff_approx)/(2*pi);