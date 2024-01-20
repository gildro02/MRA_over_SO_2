
L=12;
rho=0.4+0.2*rand(L,1);
rho=rho/sum(rho);
%{
sigma=0.3;
N_s=floor(logspace(2,6,100));
n=5;
X=10*rand(L,n);

ErrorVec_N=zeros(length(N_s),1);
for j1=1:length(N_s)
    N=N_s(j1);
    MRA_Multi_Ring
    ErrorVec_N(j1)=X_error;
end
figure
subplot(3,1,1)
loglog(N_s,ErrorVec_N)

sigma_s=logspace(-5,0.5,100);
N=1e5;
n=5;



ErrorVec_sigma=zeros(length(sigma_s),1);
for j2=1:length(sigma_s)
    sigma=sigma_s(j2);
    MRA_Multi_Ring
    ErrorVec_sigma(j2)=X_error;
end
subplot(3,1,2)
loglog(sigma_s,ErrorVec_sigma)
%}
figure
sigma=0.00;
N=1e1;
n_s=1:1:20;


ErrorVec_n=zeros(length(n_s),1);
num_rep=30;
temp_n=zeros(length(n_s),num_rep);

for j3=1:length(n_s)
    for rep=1:num_rep
    n=n_s(j3);
    X=10*rand(L,n);
    MRA_Multi_Ring
    temp_n(j3,rep)=X_error;
    end
end
ErrorVec_n=median(temp_n,2);
subplot(3,1,3)
plot(n_s,ErrorVec_n)