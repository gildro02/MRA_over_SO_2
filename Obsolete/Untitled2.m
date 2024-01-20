addpath(genpath('basis'))
addpath(genpath('common'))
addpath(genpath('examples'))
addpath(genpath('fourier'))
addpath(genpath('projections'))
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
end