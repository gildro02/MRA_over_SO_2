function [err_squared,theta_min] = circ_error_continuous_unrestricted(x,y,B)
%x,y are colomn vectors.
d_theta=1e-5;
theta=0:d_theta:(2*pi-d_theta);
[err,index]=min(vecnorm(x.*exp(1i.*(-B:B).'.*theta)-y,2,1)); %on dim 1 %monitor change to 2,1 from 2
theta_min=theta(index);
err_squared=err^2;
end
