function [x,rho] = MRA_Single_Ring(v_p,M1,M2,sigma)
L=length(M1);
W=(1/sqrt(L))*dftmtx(L);
Q=W*diag(v_p)*W';
M2_wave=Q*M2*Q';
[V,D]=eig(M2_wave,"vector");
[~,u_index]=max(min(abs(D-D.')+diag(inf(1,L)),[],2)); 
%u_index
%D
u=V(:,u_index);

x=1/sqrt(L)*ifft((fft(u)./abs(fft(u)))./v_p);
if sum(x)<0
    x=-x;
end
C_x=circul(x);
rho=C_x\M1;



end