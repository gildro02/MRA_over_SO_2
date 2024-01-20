%% Parameter Initialization
sigma=0.1;
N=1000000;
n=5;
L=12;
r=1:n;
% initialize X:
X=cell(n,1);
for m=1:n
    X{m}=(1:r(m)*L).';
end
FFT_X=cell(n,1);
for m=1:n
    FFT_X{m}=fft(X{m});
end

% initialize rho:
rho=(1:L).';
rho=rho./sum(rho);
%% Generating MRA Data
rotations=randsample(0:L-1,N,true,rho);
data=zeros(L*sum(r),N);
%data_cells=cell(n,N);
r_sum=cumsum(r);
r_start_index=[0 cumsum(r(1:end-1))];
for k=1:N
    for m=1:n
        data((L*r_start_index(m)+1):(L*r_sum(m)),k)=circshift(X{m},r(m)*rotations(k));
    end
end
noise=sigma*randn(L*r_sum(end),N);
data=data+sigma*randn(L*r_sum(end),N);
%% Generate Moments
circulant_vec=zeros(L*r_sum(end),L);
for m=1:n
    circulant_vec((L*r_start_index(m)+1):(L*r_sum(m)),:)=circul_semi(X{m},L);
end
M_1_anal=circulant_vec*rho;
M_1_emp=mean(data,2);

M_2_emp=(1/N)*(data*data.');
M_2_anal=circulant_vec*diag(rho)*circulant_vec.';

max(abs(M_1_anal-M_1_emp))
max(abs(M_2_anal-M_2_emp),[],"all")

%% M_2_wave
P_X=cell(n,1);
V=cell(n,1);
for m=1:n
    P_X{m}=abs(FFT_X{m}).^2;
    V{m}=1./sqrt(P_X{m});
end
DFT_matricies=cell(n,1);
for m=1:n
    DFT_matricies{m}=(1/sqrt(r(m)*L))*dftmtx(r(m)*L);
end
M_2_emp_wave=zeros(L*r_sum(end),L*r_sum(end));
M_2_emp_undiased=M_2_emp-(sigma^2)*eye(L*r_sum(end));
for m1=1:n
    for m2=1:n
        M_2_emp_wave((L*r_start_index(m1)+1):(L*r_sum(m1)),(L*r_start_index(m2)+1):(L*r_sum(m2)))=...
            DFT_matricies{m1}'*diag(V{m1})*DFT_matricies{m1}*...
            M_2_emp_undiased((L*r_start_index(m1)+1):(L*r_sum(m1)),(L*r_start_index(m2)+1):(L*r_sum(m2)))*...
            DFT_matricies{m2}'*diag(V{m2})*DFT_matricies{m2};
    end
end
M_2_anal_wave=zeros(L*r_sum(end),L*r_sum(end));
for m1=1:n
    for m2=1:n
        M_2_anal_wave((L*r_start_index(m1)+1):(L*r_sum(m1)),(L*r_start_index(m2)+1):(L*r_sum(m2)))=...
            DFT_matricies{m1}'*diag(V{m1})*DFT_matricies{m1}*...
            M_2_anal((L*r_start_index(m1)+1):(L*r_sum(m1)),(L*r_start_index(m2)+1):(L*r_sum(m2)))*...
            DFT_matricies{m2}'*diag(V{m2})*DFT_matricies{m2};
    end
end
M_2_anal_wave=real(M_2_anal_wave+M_2_anal_wave.')/2; %force real, symmetric
M_2_emp_wave=real(M_2_emp_wave+M_2_emp_wave.')/2; %force real, symmetric
[U_anal_top,D_anal_top]=eigs(M_2_anal_wave,L);
[U_emp_top,D_emp_top]=eigs(M_2_emp_wave,L);
D_anal_top=diag(D_anal_top);
D_emp_top=diag(D_emp_top);
[~,u_anal_index]=max(min(abs(D_anal_top-D_anal_top.')+diag(inf(L,1)),[],1),[],2);
u_anal=U_anal_top(:,u_anal_index);
[~,u_emp_index]=max(min(abs(D_emp_top-D_emp_top.')+diag(inf(L,1)),[],1),[],2);
u_emp=U_emp_top(:,u_emp_index);
X_wave_anal=cell(n,1);
for m=1:n
    X_wave_anal{m}=ifft(sqrt(P_X{m}).* ...
        fft(u_anal(L*r_start_index(m)+1:(L*r_sum(m)))) ...
        /norm(u_anal(L*r_start_index(m)+1:(L*r_sum(m)))));
end
X_wave_emp=cell(n,1);
for m=1:n
    X_wave_emp{m}=ifft(sqrt(P_X{m}).* ...
        fft(u_emp(L*r_start_index(m)+1:(L*r_sum(m)))) ...
        /norm(u_emp(L*r_start_index(m)+1:(L*r_sum(m)))));
end

