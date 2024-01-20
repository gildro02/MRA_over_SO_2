

close all
picture=imread('flower-1.jpg');
Full_Picture=im2double(rgb2gray(imresize(picture,103/1000)));
L=size(Full_Picture,1);
N=floor(L/2);
[x, y]=meshgrid(-N:N, -N:N);
r=sqrt(x.^2+y.^2);
Compressed_Picture=Full_Picture(r<=N);

%{
% Generate Basis:
[ Phi_ns, ang_freqs, rad_freqs, R_ns ]=Bessel_ns_v5( N );
Phi_ns_mat=[Phi_ns, conj(Phi_ns(:, ang_freqs~=0))];

coeff = Phi_ns_mat\Compressed_Picture;
%}
full_ang_freqs=[ang_freqs;-ang_freqs(ang_freqs~=0)];

Compressed_Proj_Picture=Phi_ns_mat*coeff;

Proj_Picture=zeros(L);
Proj_Picture(r<=N)=Compressed_Proj_Picture;


%alpha=pi;
%coeff_rot=coeff.*exp(1i*alpha.*full_ang_freqs);
%Rotated_Picture=zeros(L);
%Rotated_Picture(r<=N)=Phi_ns_mat*coeff_rot;

%Rotated_Picture=imrotate(Full_Picture,180)
%just make a LOOP%%%%%%%%%%%%%
%or mabye put all photos to cells and make a rotation cellfun


%{
figure
imagesc(real(Rotated_Picture))
title('rotated')
figure
imagesc(real(Proj_Picture))
title('projected')
%}