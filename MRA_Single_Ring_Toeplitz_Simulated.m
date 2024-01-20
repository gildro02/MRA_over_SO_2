temp=200;
a=exp(-(1/2)*sqrt((1:temp+1).'));%+1i*[0;ones(B-2,1)]; % k=0 to B
a_symm=[conj(flip(a(2:end))); a];
B=length(a)-1;
seed_hat=exp(1i.*(1:(B*2/5)).').*exp(-((1:(B*2/5)).^(1/4).')); % k=1 to B/2
%better to be more or less the same size of coeff;
[rho_hat_symm,rho_func]=Generate_Rho_Squared(seed_hat);

%rho_hat_symm=(rho_hat_symm & 1).*(1/(2*pi));
%rho_func=@(x) (1/1000)*normpdf(1000*x);

rho_hat_symm=[zeros(B-2*length(seed_hat),1);
    rho_hat_symm;
    zeros(B-2*length(seed_hat),1)];
rho_hat=rho_hat_symm((end+1)/2:end);
rho_hat_extra=zeros(B,1); % B elements
rho_hat_full=[rho_hat;rho_hat_extra];


k=(-B:B).';
W=(1/sqrt(2*B+1))*dftmtx(2*B+1);

k_half=(0:B).';
M_1_true=2*pi*a_symm.*rho_hat_symm;
M_2_true=2*pi*diag(a_symm)*toeplitz(conj(rho_hat_full))*diag(conj(a_symm));

x_time=ifft(a_symm);
rho_time=ifft(rho_hat_full);
M_1_time_true=ifft(M_1_true)/(2*pi);
M_2_time_true=W'*M_2_true*W/(2*pi);
v_p=(1/sqrt(2*B+1))./abs(fft(x_time));
[x_time_est_true,rho_time_est_true]=MRA_Single_Ring(v_p,M_1_time_true,M_2_time_true,0);


d_theta=1e-5;
theta=0:d_theta:2*pi-d_theta;
N=1e4;
sigma=0; %noise
rotations=randsample(theta,N,true,rho_func(theta)).';
a_rot=a.*exp(1i.*k_half.*rotations.')...
    +(sigma/sqrt(2*pi))*randn(B+1,N);%N colomns of rotated a vectors
a_rot_symm=[conj(flip(a_rot(2:end,:),1)) ; a_rot];

M_1_emp=mean(a_rot_symm,2);
M_2_emp=(1/N)*(a_rot_symm*a_rot_symm')-(sigma^2/(2*pi))*(k==k.'|k==-k.');
%M_2_emp=M_2_emp.*toeplitz([ones(1,B+1) zeros(1,B)]); %set zeros where we
%know there are zeros

M_1_time_emp=ifft(M_1_emp)/(2*pi);
M_2_time_emp=W*M_2_emp*W'/(2*pi);
v_p=(1/sqrt(2*B+1))./abs(fft(x_time));
[x_time_est_emp,rho_time_est_emp]=MRA_Single_Ring(v_p,M_1_time_emp,M_2_time_emp,sigma);
[err, x_time_est_emp_matched]=Circular_Error(x_time_est_emp,x_time);
err_comp=err/norm(x_time)
close all
figure
plot(0:1e-3:2*pi,rho_func(0:1e-3:2*pi))