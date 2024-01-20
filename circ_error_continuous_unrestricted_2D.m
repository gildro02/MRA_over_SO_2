% Input: x,y are column vectors of the same length Q(2B+1)x1
% Output: 1) The circular squared error between them, assuming that they represent the
% Fourier coefficients of a 2D function, defined on the unit disk. This rotation
% is unrestricted, with no quantization limit.
% 2) theta_min, the angle of rotation of the first vector where the mininum
% distance is met.

function [err_squared,theta_min] = circ_error_continuous_unrestricted_2D(x,y,Q,B)
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

