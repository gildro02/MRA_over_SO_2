function [err_squared,theta_min] = circ_error_continuous_unrestricted_2D(x,y,Q,B)
%x,y are colomn vectors of length Q(2B+1)x1
M=2*B+1;
d_theta=1e-5;
theta=0:d_theta:(2*pi-d_theta);
mat_x=mat(x,Q,B);
mat_y=mat(y,Q,B);
exp_multi_matrix=repmat(exp(1i.*(-B:B).*permute(theta,[1,3,2])),[Q,1,1]);
[err,index]=min(vecnorm(reshape(exp_multi_matrix.*mat_x-mat_y,[Q*(2*B+1),1,length(theta)]),2,1),[],3);
theta_min=theta(index);
err_squared=err^2;
end

