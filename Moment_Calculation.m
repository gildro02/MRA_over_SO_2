close all
addpath(genpath('basis'))
addpath(genpath('common'))
addpath(genpath('examples'))
addpath(genpath('fourier'))
addpath(genpath('projections'))
%% Initialyzing - comment after first run

picture=imread('flower-1.jpg');
Full_Picture=im2double(rgb2gray(imresize(picture,51/1000)));
L=size(Full_Picture,1);
N=floor(L/2);
[x, y]=meshgrid(-N:N, -N:N);
r=sqrt(x.^2+y.^2);
Compressed_Picture=Full_Picture(r<=N);
Full_Picture(r>N)=0;

% Generate Basis:
%{
[ Phi_ns, ang_freqs, rad_freqs, R_ns ]=Bessel_ns_v5( N );
Phi_ns_mat=[Phi_ns, conj(Phi_ns(:, ang_freqs~=0))];
coeff = Phi_ns_mat\Compressed_Picture;
full_ang_freqs=[ang_freqs;-ang_freqs(ang_freqs~=0)];

Compressed_Proj_Picture=Phi_ns_mat*coeff;
Proj_Picture=zeros(L);
Proj_Picture(r<=N)=Compressed_Proj_Picture;
Proj_Picture=real(Proj_Picture); %%picture we will use
%}

%% Generate Data

P=1000; %num of rotated photos;
Data=repmat(Proj_Picture,[1,1,P]);
rot_angles=360*rand(P,1); %assuming uniform distribution
%angles are in degrees


sigma=20; %noise level
for m=1:P %rotate picture by all angles.
    tic
    Data(:,:,m)=imrotate(Data(:,:,m),rot_angles(m),'bicubic','crop'); %most accurate spin
    %Data(:,:,m)=fastrotate(Data(:,:,m),rot_angles(m)); %slower, better
    Data(:,:,m)=Data(:,:,m)+sigma*randn(L).*(r<=N); %add noise
    if mod(m,100)==0
        toc
        m
        tic
    end
end
toc

"coeff of data"
%calculate noise:
coeff_Noise_positive=sigma*randn(length(ang_freqs),P);
coeff_Noise=[coeff_Noise_positive; coeff_Noise_positive(ang_freqs~=0,:)];
%analitic coefficients (rotated pictures) + real noise:
coeff_Data=repmat(coeff,1,P).*exp(1i*full_ang_freqs*rot_angles'*2*pi/360) + coeff_Noise;
coeff_mean_analitic=cumsum(coeff_Data,2)./repmat(1:P,[length(coeff) 1]);
%%
%calculate cascading mean;
MeansVec=bsxfun(@rdivide,cumsum(Data,3),permute(1:P,[1 3 2]));

%{
MeansVec=zeros(L,L,P);
MeansVec(:,:,1)=Data(:,:,1);
for m=2:P
    tic
    MeansVec(:,:,m)=((m-1)*MeansVec(:,:,m-1)+Data(:,:,m))/m;
    if mod(m,100)==0
        toc
        m
        tic
    end
end
toc
%}

Compressed_MeansVec=reshape(MeansVec(repmat(r<=N,[1,1,P])),[length(Compressed_Picture),P]);
tic
"big process"
coeff_mean = Phi_ns_mat\Compressed_MeansVec;
toc

Zero_Freq_Mean=coeff_mean_analitic(full_ang_freqs==0,:);
Zero_Freq_Original=coeff(full_ang_freqs==0,:);
ErrorVec=(vecnorm((Zero_Freq_Mean-Zero_Freq_Original)./Zero_Freq_Original)); %add back sqrt(1:P)
figure
loglog(1:P,ErrorVec); %Error line
hold on
loglog(1:P,(1:P).^(-1/2)); %slope 1/2 line







bound=P/2:P;
boundVec=log(bound);
newErrorVec=log(ErrorVec(bound));
%need to see -1/2 slope line!!

%{
%brute force rotate
alpha=pi;
coeff_rot=coeff.*exp(1i*alpha.*full_ang_freqs);
Rotated_Picture=zeros(L);
Rotated_Picture(r<=N)=Phi_ns_mat*coeff_rot;
%}