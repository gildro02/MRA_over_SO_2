%{


picture=imread('flower-1.jpg');
pic_size=1000;
Num_picture=im2double(rgb2gray(imresize(picture,pic_size/size(picture,1))));

imagesc(Num_picture)

R=floor(pic_size/2);
[x, y]=meshgrid(-R:R, -R:R);
r=sqrt(x.^2+y.^2);
%}
clear all
c=30; %group size 
m=2; %degeneration level
B=10; %bandwidth
%c must be >=2B+1
if ~(c>=2*B+1)
    error("c must be >=2B+1")
end
N=500; %number of pictures
sigma=0; %noise level
omega_s=m*c; %sampling rate=amount of samples on [0,2pi]
T_s=2*pi/omega_s; %sampling period

x=@(t) 2*pi-t; %between 0 and 2pi
[x_proj,coeff]=Project_Onto_Exponents(x,B);
%x=@(t) sin(t)-2*cos(2*t);
rho=1:c;
rho=rho./sum(rho);
x_sample=real(x_proj(0:T_s:2*pi-T_s));

data=Generate_MRA_Data(rho,N,sigma,x_sample,0:m:m*(c-1));
data_by_remainder=zeros(m,c,N);
for j=1:N
    data_by_remainder(:,:,j)=reshape(data(:,j),[m,c]);
end

M1_by_remainder=mean(data_by_remainder,3).';
data_by_remainder_permuted=permute(data_by_remainder,[2,3,1]);
M2_by_remainder=zeros(c,c,m);
for j=1:m
    M2_by_remainder(:,:,j)=data_by_remainder_permuted(:,:,j)*data_by_remainder_permuted(:,:,j).'/N-sigma^2*eye(c);
end

v_p_by_remainder=(1/sqrt(c))./abs( ...
    fft( ... %fft by colomn
    reshape(x_sample,[m,c]).')); %reshape to seperate to different remainders, transpose to make those colomn vectors


x_est_by_remainder=zeros(c,m);
rho_est_by_remainder=zeros(c,m);
for j=1:m
    [x_est_by_remainder(:,j),rho_est_by_remainder(:,j)]=MRA_Single_Ring(v_p_by_remainder(:,j),M1_by_remainder(:,j),M2_by_remainder(:,:,j),sigma);
end
x_est_by_remainder=real(x_est_by_remainder);
%v=zeros(c,2*B+1,m);
v=exp(1i.*(-B:B).*...
    ((2*pi/(c*m)).*permute(0:m-1,[1,3,2])+...
    (2*pi/c).*(0:c-1).'));

alpha_est_by_remainder=zeros(2*B+1,m);
for j=1:m
    alpha_est_by_remainder(:,j)=pinv(v(:,:,j))*x_est_by_remainder(:,j);
end
phase=alpha_est_by_remainder(B+1:end,:)./coeff
log_phase=log(alpha_est_by_remainder(B+1:end,:)./coeff)
figure
plot(0:B,mod(real(log_phase/1i),2*pi))
circshift(phase,1,1)
%need to decompose to different remainders
%positive_coeff=