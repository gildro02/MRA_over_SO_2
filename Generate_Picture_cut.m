% Input: 1) name: a char vector of the the file name of the picture to be
% spanned by Fourier-Bessel Coefficient.
% 2) B,Q: the requested cut-off frequency for the angular and radial
% frequencies. all freqiencies above these will be thrown out.
% 3) isUniformPowerSpectrum: a binary value. If true, the function ensures
% that the power spectrum is uniform, by normalizing each F-B coefficient
% to be of size 1. If false, it does not alter the power spectrum after
% cutting excessive frequencies.
% 
% Output: 1) coeff: the F-B coefficient vector, organized first by angular
% and then by radial frequencies (as a dictionary).
% 2) Proj_Picture: the picture corresponding to those F-B coefficients.
% 3) Phi_ns_mat: a matrix with its column vectors being the basis vectors
% of the F-B basis, orginized in the same way as coeff.

function [coeff,Proj_Picture,Phi_ns_mat]=Generate_Picture_cut(name,B,Q,image_size,isUniformPowerSpectrum)
%close all
addpath(genpath('basis'))
addpath(genpath('common'))
addpath(genpath('projections'))

%if ~exist('picture','var')
picture=imread(name);
Full_Picture=im2double(rgb2gray(imresize(picture,image_size./size(picture,1))));
L=size(Full_Picture,1);
N=floor(L/2);
[x, y]=meshgrid(-N:N, -N:N);
r=sqrt(x.^2+y.^2);
Compressed_Picture=Full_Picture(r<=N);
Full_Picture(r>N)=0;

% Generate Basis:
[ Phi_ns, ang_freqs, rad_freqs, R_ns ]=Bessel_ns_v5( N );
%end
k_symm_1B=-B:B;
q_vec=1:Q; %odd, no start from 0
ang_freqs_cut=repelem(k_symm_1B,Q);
rad_freqs_cut=reshape(repmat(q_vec.',[1,2*B+1]),[1,Q*(2*B+1)]);

Phi_ns_mat=zeros(length(Phi_ns(:,1)),Q*(2*B+1));
for k_index=-B:B
    for q_index=1:Q
        if k_index>=0
            Phi_ns_mat(:,ang_freqs_cut==k_index & rad_freqs_cut==q_index)...
                = Phi_ns(:,ang_freqs==k_index & rad_freqs==q_index);
        else
            Phi_ns_mat(:,ang_freqs_cut==k_index & rad_freqs_cut==q_index)...
            = conj(Phi_ns(:,ang_freqs==-k_index & rad_freqs==q_index));
        end
    end
end

coeff = Phi_ns_mat\Compressed_Picture;
if isUniformPowerSpectrum
    coeff=coeff./abs(coeff); %ensure uniform power spectrum
end
Compressed_Proj_Picture=Phi_ns_mat*coeff;
Proj_Picture=zeros(L,L);
Proj_Picture(r<=N)=Compressed_Proj_Picture;
Proj_Picture=real(Proj_Picture); %%picture we will use


%figure
%imagesc(Proj_Picture)
%colorbar