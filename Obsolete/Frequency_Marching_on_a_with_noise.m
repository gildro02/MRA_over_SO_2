temp=20;
a=exp(-(1/2)*sqrt((1:temp-1).'));%+1i*[0;ones(B-2,1)]; % k=0 to B
a_symm=[conj(flip(a(2:end))); a];
B=length(a)-1;
seed_hat=exp(-(1:B/2).'); % k=1 to B/2
[rho_hat_symm,rho_func]=Generate_Rho_Squared(seed_hat);
rho_hat=rho_hat_symm((end+1)/2:end);
rho_hat_extra=zeros(B,1); % B elements
rho_hat_full=[rho_hat;rho_hat_extra];

k=(-B:B).';
k_half=(0:B).';
M_1_true=2*pi*a_symm.*rho_hat_symm;
M_2_true=2*pi*diag(a_symm)*toeplitz(conj(rho_hat_full))*diag(conj(a_symm));

d_theta=1e-5;
theta=0:d_theta:2*pi-d_theta;
N=1e6;
sigma=6; %noise
rotations=randsample(theta,N,true,rho_func(theta)).';
a_rot=a.*exp(1i.*k_half.*rotations.')...
    +(sigma/sqrt(2*pi))*randn(B+1,N);%N colomns of rotated a vectors
a_rot_symm=[conj(flip(a_rot(2:end,:),1)) ; a_rot];

M_1_emp=mean(a_rot_symm,2);
M_2_emp=(1/N)*(a_rot_symm*a_rot_symm')-(sigma^2/(2*pi))*(k==k.'|k==-k.');
%norm(M_2_emp-M_2_true,"fro")
%norm(M_1_emp-M_1_true,"fro")

M_2_wave=M_2_emp.*toeplitz([1./conj(M_1_emp(k>=0)) ; zeros(B,1)]);

weights=cell(1,B-1); %max number of approximations x max freq-1
for n=1:B-1
    weights{n}=exp(-abs(ceil(n/2)-(1:n).')); %lorentzian of length n centered at ceil(n/2)
    weights{n}=weights{n}./sum(weights{n}); %normalize to sum 1.
end


a_approx=zeros(B+1,1);
a_approx(k_half==0)=M_1_emp(k==0);
a_approx(k_half==1)=sqrt(M_2_wave(k==1,k==1)*M_1_emp(k==0));
a_approx_abs=sqrt(a_approx(k_half==0).*diag(M_2_wave(k>=0,k>=0)));

for n=2:B
    a_approx_raw=M_2_wave(k==n,(1<=k)&(k<n)).'.*...
                 flip(a_approx((1<=k_half)&(k_half<n)))./...
                 conj(a_approx((1<=k_half)&(k_half<n)));
    %regular mean
    %mean_a_approx_raw=mean(a_approx_raw);

    %custom mean:
    mean_a_approx_raw=sum(a_approx_raw.*weights{n-1});

    %using known absolute value from M2_wave:
    a_approx_phase=mean_a_approx_raw./abs(mean_a_approx_raw);
    a_approx(k_half==n)=a_approx_phase*a_approx_abs(k_half==n);
    %not using known absolute value from M2_wave:
    %a_approx(k_half==n)=mean_a_approx_raw;
end