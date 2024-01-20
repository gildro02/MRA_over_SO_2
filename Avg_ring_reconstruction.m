L=10; %length of signal
R=20; %num of rings (different signals)
N=500; %num of noisy copies of the same signal.
sigma=0;
x_s=100*rand(L,R)+sigma*randn(L,R); %colomns are different rings (signals);
%x_s=repmat((1:L).',[1,R]);
rho=rand(L,1);
%rho=(1:L).';
rho=rho/sum(rho);

M1_true=ifft(repmat(fft(rho),[1,R]).*fft(x_s));
M2_true=zeros(L,L,R);
for r=1:R
M2_true(:,:,r)=circul(x_s(:,r))*diag(rho)*circul(x_s(:,r)).';
end

v_p=(1/sqrt(L))./abs(fft(x_s));
%v_p=1./abs(fft(x_s));

data=zeros(L,N,R);
for r=1:R
data(:,:,r)=Generate_MRA_Data(rho,N,sigma,x_s(:,r));
end

M1_emp=permute(mean(data,2),[1,3,2]);
M2_emp=zeros(L,L,R);
for r=1:R
M2_emp(:,:,r)=data(:,:,r)*(data(:,:,r)).'/N -(sigma^2)*eye(L);
end




x_est=zeros(L,R);
rho_est=zeros(L,R);
error=zeros(1,R);
for r=1:R
[x_est(:,r),rho_est(:,r)]=MRA_Single_Ring(v_p(:,r),M1_true(:,r),M2_true(:,:,r),sigma);
error(r)=Circular_Error(x_est(:,r),x_s(:,r));
end
max_error=max(error)