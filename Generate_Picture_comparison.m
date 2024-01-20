%close all
addpath(genpath('basis'))
addpath(genpath('common'))
addpath(genpath('examples'))
addpath(genpath('fourier'))
addpath(genpath('projections'))

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
k_full=[-flip(ang_freqs(ang_freqs~=0)); ang_freqs];
q_full=
Phi_ns_mat=[flip(conj(Phi_ns(:, ang_freqs~=0)),2), Phi_ns];
full_ang_freqs=[-flip(ang_freqs(ang_freqs~=0)); ang_freqs];
full_rad_freqs=[flip(rad_freqs(ang_freqs~=0)); rad_freqs];
%this approach also flips the q's which is unintentional.

ang_freqs_max=max(ang_freqs);
B=10;
Q_vec=zeros(B+1,1); %number of radial frequencies per angular frequency
                                %from 0 to B
for index1=0:ang_freqs_max
    Q_vec(index1+1)=sum(ang_freqs==index1);
end
end
coeff = Phi_ns_mat\Compressed_Picture;

%set overflow coefficients to be zero, so that Q_k is constant
B=20;
coeff(abs(full_ang_freqs)>=B|full_rad_freqs>=max(full_rad_freqs(full_ang_freqs==B)))=0;

Compressed_Proj_Picture=Phi_ns_mat*coeff;
Proj_Picture=zeros(L);
Proj_Picture(r<=N)=Compressed_Proj_Picture;
Proj_Picture=real(Proj_Picture); %%picture we will use

toc

figure
imagesc(Proj_Picture)