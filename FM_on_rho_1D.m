sigma_vec=logspace(-2,2,50);
err_squared_a=zeros(length(sigma_vec),1);
B=10;
N=1e6;

%a=exp(2*pi*1i*rand(B,1));%.*exp(-(1/2)*sqrt((1:B).'));%+1i*[0;ones(B-2,1)]; % k=0 to B
a=[0.707671873258695 + 0.706541237153593i;0.988238403497549 + 0.152921083741304i;-0.249817202634132 + 0.968293016223941i;-0.411632223458529 + 0.911350049437968i;-0.568871246227187 - 0.822426595639955i;0.963616145197972 - 0.267289963735642i;0.919568728484166 - 0.392929196667815i;-0.965195103437267 + 0.261530901234871i;0.0597903158092206 + 0.998210958733390i;0.0872124228383627 - 0.996189737601559i];
a_symm_1B=[conj(flip(a(1:end))); 1 ;a];
a_half_1B= [1;a];
k_symm_1B=-B:B;
k_half_1B=0:B;
k_symm_2B=-2*B:2*B;
k_half_2B=0:2*B;

temp=[0.00195583374369666 + 0.00530869159003379i;0.00195583374369666 + 0.00474988194897760i;0.00195583374369666 + 0.00419107230792141i;0.00195583374369666 + 0.00363226266686522i;0.00195583374369666 + 0.00307345302580903i;0.00195583374369666 + 0.00251464338475285i;0.00195583374369666 + 0.00195583374369666i;0.00195583374369666 + 0.00139702410264047i;0.00195583374369666 + 0.000838214461584282i;0.00195583374369666 + 0.000279404820528094i;0.00195583374369666 - 0.000279404820528094i;0.00195583374369666 - 0.000838214461584282i;0.00195583374369666 - 0.00139702410264047i;0.00195583374369666 - 0.00195583374369666i;0.00195583374369666 - 0.00251464338475285i;0.00195583374369666 - 0.00307345302580903i;0.00195583374369666 - 0.00363226266686522i;0.00195583374369666 - 0.00419107230792141i;0.00195583374369666 - 0.00474988194897760i;0.00195583374369666 - 0.00530869159003379i]...
    .*exp(1i*(1*rand(2*B,1)).^0.2)*3;
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

weights=cell(1,B); %max number of approximations x max freq-1
for k=2:B
    weights{k}=exp(-abs(ceil(k/2)-(1:k-1).')); %lorentzian of length n centered at ceil(n/2)
    weights{k}=weights{k}./sum(weights{k}); %normalize to sum 1.
end

for index=1:length(sigma_vec)
sigma=sigma_vec(index); %noise

M_1_true=2*pi*a_symm_1B.*rho_hat_symm_1B;
T=toeplitz(conj(rho_hat_half_2B));
M_2_true=2*pi*diag(a_symm_1B)*T*diag(conj(a_symm_1B))+...
    (sigma^2/(2*pi))*(k_symm_1B==k_symm_1B.'|k_symm_1B==-k_symm_1B.');

rotations=randsample(t,N,true,abs(rho_func(t))).';
%absolute value to avoid numeric negatives, rho should always be positive.
a_rot=a_half_1B.*exp(1i.*k_half_1B.'.*rotations.')...
    +(sigma/sqrt(2*pi))*randn(B+1,N);%N colomns of rotated a vectors
a_rot_symm=[conj(flip(a_rot(2:end,:),1)) ; a_rot];

M_1_emp=mean(a_rot_symm,2);
M_2_emp=(1/N)*(a_rot_symm*a_rot_symm');

moment_type="emp";
%moment_type="true";
if strcmp(moment_type,"emp")==1
    M_1=M_1_emp;
    M_2=M_2_emp;
end
if strcmp(moment_type,"true")==1
    M_1=M_1_true;
    M_2=M_2_true;
end
M_2_wave = M_2-(sigma^2/(2*pi))*(k_symm_1B==k_symm_1B.'|k_symm_1B==-k_symm_1B.');

S=2*pi*diag(1./M_1)*M_2_wave*diag(1./M_1)';

rho_approx_half_1B=zeros(length(k_half_1B),1);
rho_approx_tilde_half_1B= zeros(length(k_half_1B),1);

rho_abs_half_1B = sqrt(abs(1./(2*pi*diag(S(k_symm_1B>=0,k_symm_1B>=0)))));
rho_approx_half_1B(k_half_1B==0)=1/(2*pi);
%arbitrarily pick phase of 1'st coeff to be 1:
rho_approx_half_1B(k_half_1B==1) = rho_abs_half_1B(k_half_1B==1); 

for k=2:B
   rho_approx_tilde_half_1B(k_half_1B==k) = sum(weights{k}.*flip(rho_approx_half_1B(1<=k_half_1B & k_half_1B <= k-1))./...
       conj(rho_approx_half_1B(1<=k_half_1B & k_half_1B <= k-1))./...
       S(k_symm_1B==k, 1<=k_symm_1B & k_symm_1B <= k-1).');
   
   rho_approx_half_1B(k_half_1B==k)=rho_abs_half_1B(k_half_1B==k).*...
       rho_approx_tilde_half_1B(k_half_1B==k)./abs(rho_approx_tilde_half_1B(k_half_1B==k));
end
rho_approx_symm_1B= [conj(flip(rho_approx_half_1B(2:end))); rho_approx_half_1B];
%[err_squared_rho,theta_min]=circ_error_continuous_unrestricted(rho_approx_symm_1B,rho_hat_symm_1B,B)

a_approx_symm_1B=M_1./(2*pi.*rho_approx_symm_1B);
[err_squared_a(index),theta_min]=circ_error_continuous_unrestricted(a_approx_symm_1B,a_symm_1B,B);

end
figure
loglog(sigma_vec,err_squared_a)