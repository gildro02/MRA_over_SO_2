close all
addpath(genpath('basis'))
addpath(genpath('common'))
addpath(genpath('examples'))
addpath(genpath('fourier'))
addpath(genpath('projections'))
%% Initialyzing - comment after first run
tic
if ~exist('picture','var')
picture=imread([pwd '/Flower_Images' '/flower-1.jpg']);
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
toc
%% Generate Data
B=20; %Band Limit
n=B/2;
seed_vec=1i./(1:B/2).^2;
[Distribution_coeff,Distribution]=Generate_Rho_Squared(seed_vec);
Distribution_coeff_indecies=-B:B;
P=1000; %num of rotated photos;
Data=repmat(Proj_Picture,[1,1,P]);
dx=1e-5; %Quantization of [0,2pi].
rot_angles=randsample(0:dx:2*pi,P,true,dx*Distribution(0:dx:2*pi))'; %randomize from distribution[0,2pi].

sigma=1; %noise level
%{
for m=1:P %rotate picture by all angles.
    tic
    Data(:,:,m)=imrotate(Data(:,:,m),(360/(2*pi))*rot_angles(m),'bicubic','crop'); %most accurate spin
    %Data(:,:,m)=fastrotate(Data(:,:,m),(360/(2*pi))*rot_angles(m)); %slower, better
    Data(:,:,m)=Data(:,:,m)+sigma*randn(L).*(r<=N); %add noise
    if mod(m,100)==0
        toc
        m
        tic
    end
end
toc
%}

"coeff of data"
%calculate noise:
coeff_Noise_positive=sigma*randn(length(ang_freqs),P);
coeff_Noise=[coeff_Noise_positive; coeff_Noise_positive(ang_freqs~=0,:)];
%analitic coefficients (rotated pictures) + real noise:
coeff_Data=repmat(coeff,1,P).*exp(-1i*full_ang_freqs*(rot_angles')) + coeff_Noise;

%% Calculate Moments (Mean + Auto Corrolation)
%legacy operations in space domain
%{
%calculate cascading mean;
MeansVec=bsxfun(@rdivide,cumsum(Data,3),permute(1:P,[1 3 2]));


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


Compressed_MeansVec=reshape(MeansVec(repmat(r<=N,[1,1,P])),[length(Compressed_Picture),P]);
tic
"big process"
coeff_mean = Phi_ns_mat\Compressed_MeansVec;
toc
%}
%legacy zero-freqs for uniform distribution
%{
%zero frequencies:
Zero_Freq_Mean=coeff_mean_analitic(full_ang_freqs==0,:);
Zero_Freq_Original=coeff(full_ang_freqs==0);
Error_Vec=vecnorm((Zero_Freq_Mean-repmat(Zero_Freq_Original,[1,P]))./Zero_Freq_Original);
%}
%all frequencies:
ang_freqs_max=max(ang_freqs);
altered_Distribution_coeff=zeros(length(coeff),1);
for index=-ang_freqs_max:ang_freqs_max
    temp=Distribution_coeff(Distribution_coeff_indecies==index);
    if ~isempty(temp)
    altered_Distribution_coeff(full_ang_freqs==index)=temp;
    end
end
%first moment
coeff_mean_Expected=coeff.*altered_Distribution_coeff.*(2*pi); %*(2pi)
coeff_mean_analitic=cumsum(coeff_Data,2)./repmat(1:P,[length(coeff) 1]);
Full_Error_Mean=vecnorm(coeff_mean_analitic-repmat(coeff_mean_Expected,[1,P]))/norm(coeff_mean_Expected);

%second moment
ang_freq_index_Mat=full_ang_freqs-full_ang_freqs.'; %colomn minus row
rad_freq_index_Mat=full_rad_freqs-full_rad_freqs.';
Distribution_coeff_Mat=zeros(length(coeff));
%place correct values to Distribution_coeff_Mat
for index=-ang_freqs_max:ang_freqs_max
    temp=Distribution_coeff(Distribution_coeff_indecies==index);
    if ~isempty(temp)
    Distribution_coeff_Mat(ang_freq_index_Mat==index)=temp;
    end
end
%Bias occurs only for q1=q2,k1=+-k2
Logical_Bias_Mat=zeros(length(coeff));
Logical_Bias_Mat((rad_freq_index_Mat==0)&ang_freq_index_Mat==0)=1;
Logical_Bias_Mat((rad_freq_index_Mat==0)&(full_ang_freqs+full_ang_freqs.'==0))=1;
coeff_autocorrelation_Expected=coeff*coeff'.*Distribution_coeff_Mat.*(2*pi)+(sigma^2)*Logical_Bias_Mat; %'=complex conjugate transpose
Full_Error_Autocorrolation=zeros(1,P);

"analitical autocorrolation calc"
coeff_autocorrelation_analitic=zeros(length(coeff)); %initialize
norm_autocorrolation=norm(coeff_autocorrelation_Expected,"fro");
time=tic;
tic
for m=1:P
    %running mean: for m=1 we get coeff_Data(:,1)*coeff_Data(:,1)';
    coeff_autocorrelation_analitic=(coeff_autocorrelation_analitic*(m-1)+coeff_Data(:,m)*coeff_Data(:,m)')/m;
    Full_Error_Autocorrolation(m)=norm(coeff_autocorrelation_analitic-coeff_autocorrelation_Expected,"fro")...
        ./norm_autocorrolation;
    if mod(m,100)==0
        toc
        m
        tic
    end
end
toc
toc(time)
%{
coeff_autocorrelation_per_iteration=zeros(length(coeff),length(coeff),P);

"coeff_autocorrelation_per_iteration"
tic
for m=1:P
    coeff_autocorrelation_per_iteration(:,:,m)=coeff_Data(:,m)*(coeff_Data(:,m)'); %complex conjugate
end
toc

%}


%% Plots
figure
loglog(1:P,Full_Error_Mean); %Error line for mean
hold on
loglog(1:P,Full_Error_Autocorrolation) %Error line for autocorrolation
hold on
loglog(1:P,(1:P).^(-1/2)); %slope -1/2 line
legend("Error of Mean","Error of Autocorrolation","Slope -1/2");

figure
plot(-2*pi:dx:2*pi,Distribution(-2*pi:dx:2*pi))





bound=P/2:P;
boundVec=log(bound);
newErrorVec=log(Full_Error_Mean(bound));
%need to see -1/2 slope line!!

%{
%brute force rotate
alpha=pi/3;
coeff_rot=coeff.*exp(-1i*alpha.*full_ang_freqs);
Rotated_Picture=zeros(L);
Rotated_Picture(r<=N)=real(Phi_ns_mat*coeff_rot);
figure
imagesc(Rotated_Picture);
%}