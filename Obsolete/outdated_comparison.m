B=10;
%a=exp(2*pi*1i*rand(B,1));%.*exp(-(1/2)*sqrt((1:B).'));%+1i*[0;ones(B-2,1)]; % k=0 to B
a=[0.707671873258695 + 0.706541237153593i;0.988238403497549 + 0.152921083741304i;-0.249817202634132 + 0.968293016223941i;-0.411632223458529 + 0.911350049437968i;-0.568871246227187 - 0.822426595639955i;0.963616145197972 - 0.267289963735642i;0.919568728484166 - 0.392929196667815i;-0.965195103437267 + 0.261530901234871i;0.0597903158092206 + 0.998210958733390i;0.0872124228383627 - 0.996189737601559i];
a_symm_1B=[conj(flip(a(1:end))); 1 ;a];

k_symm_1B=-B:B;
k_half_1B=0:B;
k_symm_2B=-2*B:2*B;
k_half_2B=0:2*B;

f=logspace(log10(500*B),log10(10*B),40).';
%f=logspace(-4,-2,40);
bound=zeros(length(f),1);
err_squared_good=zeros(length(f),1);
err_squared_bad=zeros(length(f),1);
dist_from_circ=zeros(length(f),1);
delta_kappa=zeros(length(f),1);
for index=1:length(f)
%temp=exp(1i*(1:2*B).'./(100*B));
%temp=ones(2*B,1);
%temp= exp(-(abs(1:2*B)).^(1/50).')*0+ones(2*B,1)+0.5*i*exp(-(abs(1:2*B)).^(0.01).');
%temp = ones(2*B,1) + 3i.*ones(2*B,1);
temp=[0.00195583374369666 + 0.00530869159003379i;0.00195583374369666 + 0.00474988194897760i;0.00195583374369666 + 0.00419107230792141i;0.00195583374369666 + 0.00363226266686522i;0.00195583374369666 + 0.00307345302580903i;0.00195583374369666 + 0.00251464338475285i;0.00195583374369666 + 0.00195583374369666i;0.00195583374369666 + 0.00139702410264047i;0.00195583374369666 + 0.000838214461584282i;0.00195583374369666 + 0.000279404820528094i;0.00195583374369666 - 0.000279404820528094i;0.00195583374369666 - 0.000838214461584282i;0.00195583374369666 - 0.00139702410264047i;0.00195583374369666 - 0.00195583374369666i;0.00195583374369666 - 0.00251464338475285i;0.00195583374369666 - 0.00307345302580903i;0.00195583374369666 - 0.00363226266686522i;0.00195583374369666 - 0.00419107230792141i;0.00195583374369666 - 0.00474988194897760i;0.00195583374369666 - 0.00530869159003379i]...
    .*exp(1i*(rand(2*B,1)).^0.2./f(index))*3;
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

%simulation:
d_theta=1e-5;
theta=0:d_theta:2*pi-d_theta;
N=1e5;
sigma=0.1; %noise
rotations=randsample(theta,N,true,rho_func(theta)).';
a_rot=a.*exp(1i.*k_half.*rotations.')...
    +(sigma/sqrt(2*pi))*randn(B+1,N);%N colomns of rotated a vectors
a_rot_symm_1B=[conj(flip(a_rot(2:end,:),1)) ; a_rot];

M_1_emp=mean(a_rot_symm_1B,2);
M_2_emp=(1/N)*(a_rot_symm_1B*a_rot_symm_1B')-...
    (sigma^2/(2*pi))*(k_symm_1B==k_symm_1B.'|k_symm_1B==-k_symm_1B.');

W=(1/sqrt(2*B+1))*dftmtx(2*B+1);
h=(1/(2*B+1))*(k_half_2B.'.*[0 ; conj(flip(rho_hat_half_2B(2:end)))]...
    +(2*B+1-k_half_2B.').*rho_hat_half_2B);
C_h=circul(h);
T=toeplitz(conj(rho_hat_half_2B));
M_1=2*pi*a_symm_1B.*rho_hat_symm_1B;
M_2_T=2*pi*diag(a_symm_1B)*T*diag(conj(a_symm_1B));
M_2_C=2*pi*diag(a_symm_1B)*C_h*diag(conj(a_symm_1B));
P_a_symm_1B=abs(a_symm_1B).^2;
a_symm_1B_wave=a_symm_1B./sqrt(P_a_symm_1B);

M_2_T_mathcal=(1/(2*pi))*W'*diag(1./sqrt(P_a_symm_1B))*M_2_T*diag(1./sqrt(P_a_symm_1B))*W;
M_2_C_mathcal=(1/(2*pi))*W'*diag(1./sqrt(P_a_symm_1B))*M_2_C*diag(1./sqrt(P_a_symm_1B))*W;

%norm(M_2_T_mathcal-circul(ifft(a_symm_1B_wave))*W'*T*W*circul(ifft(a_symm_1B_wave))')

M_2_T_mathcal=(M_2_T_mathcal'+M_2_T_mathcal)/2; %force hermitian
M_2_C_mathcal=(M_2_C_mathcal'+M_2_C_mathcal)/2;

[V_T,D_T]=eig(M_2_T_mathcal,"vector");
[V_C,D_C]=eig(M_2_C_mathcal,"vector");

[D_T_sorted,I_T_sort]=sort(D_T,"descend");
[D_C_sorted,I_C_sort]=sort(D_C,"descend");

V_T_sorted=V_T(:,I_T_sort);
V_C_sorted=V_C(:,I_C_sort);

[~,kappa]=max(min(abs(D_T_sorted-D_T_sorted.')+diag(inf(1,2*B+1)),[],2));
%kappa=21;
delta_kappa(index)=max(min(abs(D_T_sorted(kappa)-D_C_sorted(1:end~=kappa))),...
    min(abs(D_C_sorted(kappa)-D_T_sorted(1:end~=kappa))));

u=fft(V_T_sorted(:,kappa));
u_renormilized_bad=exp(1i* (angle(a_symm_1B(k_symm_1B==0))-angle(u(k_symm_1B==0)))) *...
    sqrt(P_a_symm_1B).* u; %./abs(u);
[err_squared_bad(index),l_bad]=circ_error_continuous(u_renormilized_bad,a_symm_1B,B);
u_renormilized_good=exp(1i* (angle(a_symm_1B(k_symm_1B==0))-angle(u(k_symm_1B==0)))) *...
    sqrt(P_a_symm_1B).* u./abs(u);
[err_squared_good(index),l_good]=circ_error_continuous(u_renormilized_good,a_symm_1B,B);

dist_from_circ(index)=sum(abs((conj(rho_hat_half_2B(2:end))-flip(rho_hat_half_2B(2:end))).^2)...
    .*(k_half_2B(2:end).'.*(flip(k_half_2B(2:end).')))./(2*B+1));

bound_rho=(2*B+1)*2*(1-sqrt(1-dist_from_circ(index)/(delta_kappa(index)^2)));
bound(index)=bound_rho*max(abs(a_symm_1B));
if imag(bound(index))~=0
    bound(index) = -1;
end

end