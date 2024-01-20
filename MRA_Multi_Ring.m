%% Parameter Initialization
%sigma=0;
%N=1000;
%n=5;
%L=12;
W=(1/sqrt(L))*dftmtx(L);

%X=10*rand(L,n);
%X=repmat((1:L).',[1 n]);

%rho=rand(L,1);
%rho=(1:L).';

%rho=rho/sum(rho);


%% Generating MRA Data

rotation=randsample(0:L-1,N,true,rho);
data=zeros(L,n,N);
for m=1:N
    data(:,:,m)=circshift(X,rotation(m),1);
end
%}
data=data+sigma*randn(L,n,N); %add noise
%% Generating Moments
M1_anal=circul(rho)*X;
M1_emp=mean(data,3);

F_X=fft(X);
v_p=1./abs(F_X);
F_X_wave=F_X.*v_p;
X_wave=ifft(F_X_wave);

circ_mat=zeros(n*L,L);
for k=1:n
    circ_mat(1+(k-1)*L:k*L,:)=circul(X_wave(:,k));
end

M2_anal=circ_mat*diag(rho)*circ_mat.';
M2_anal=(M2_anal+M2_anal')./2;
data_linear=reshape(data,L*n,N);
M2_emp_non_normilized=(1/N)*(data_linear*data_linear.');
%{
normalization_matrix=zeros(n*L,n*L);
for k=1:n
   normalization_matrix(1+(k-1)*L:k*L,1+(k-1)*L:k*L)=W*diag(v_p(:,k))*W';
end
M2_emp=normalization_matrix*(M2_emp_non_normilized-sigma^2*eye(n*L))*normalization_matrix;
%}
%equivalent:

M2_emp=zeros(n*L,n*L);
for k1=1:n
    for k2=1:n
        M2_emp(1+(k1-1)*L:k1*L,1+(k2-1)*L:k2*L)=W*diag(v_p(:,k1))*W'*...
            (M2_emp_non_normilized(1+(k1-1)*L:k1*L,1+(k2-1)*L:k2*L)-(k1==k2)*(sigma^2)*eye(L))*...
            W*diag(v_p(:,k2))*W';
    end
end
M2_emp=(M2_emp+M2_emp')/2;% force hermitian
%% Finding X and Rho
[V_top,D_top]=eigs(M2_emp,L,"largestabs");
D_top=diag(D_top);
%Find the index of the eigenvalue furthest away from all other eigenvalues
[~,u_index]=max(min(abs(D_top-D_top.')+diag(inf(L,1)),[],1),[],2);
%mabye it could help to just take the top half eigenvalues, and hope that
%one of them is non degenerate?
u=V_top(:,u_index);
%need to normalize
U_by_colomn=reshape(u,[L,n]);
F_U_by_colomn=fft(U_by_colomn);
%Renormilize in fourier domain:
F_X_wave_estimate=F_U_by_colomn./abs(F_U_by_colomn);
F_X_estimate=F_X_wave_estimate./v_p;

X_estimate=real(ifft(F_X_estimate));
%temporary sign designation:
if sum(X_estimate)<0
    X_estimate=-X_estimate;
end
X_error=Circular_Error(X_estimate,X)/n

D_emp_sorted=sort(eig(M2_emp,"vector"),"descend");
%D_emp=flip(eig(M2_emp));
%D_anal=flip(eig(M2_anal));
%D_true=flip(sort(rho*n));


%if rho has a value thats too low,the corresponding eigenvalue can sometimes
%be sucked into the 0-eigenvalues range. this is more prevalent for high
%sigma/low N. This might be solvable if we than pick an eigenvalue thats
%not one of the smaller ones, though we might encounter a problem with its
%degenerocity.

%circ_mat_X_estimate=zeros(n*L,L);
rho_estimate_mat=zeros(L,n);
for k=1:n
    %circ_mat_X_estimate(1+(k-1)*L:k*L,:)=circul(X_estimate(:,k));
    rho_estimate_mat(:,k)=circul(X_estimate(:,k))\M1_emp(:,k);
end
