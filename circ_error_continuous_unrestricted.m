% Input: x,y are column vectors of the same length 2B+1
% Output: 1) The circular squared error between them, assuming that they represent the
% Fourier coefficients of a 1D function, defined on [0,2pi]. This rotation
% is unrestricted, with no quantization limit.
% 2) theta_min, the angle of rotation of the first vector where the mininum
% distance is met.
function [err_squared,theta_min] = circ_error_continuous_unrestricted(x,y,B)
d_theta=1e-5;
theta=0:d_theta:(2*pi-d_theta);
[err,index]=min(vecnorm(x.*exp(1i.*(-B:B).'.*theta)-y,2,1)); %on dim 1 %monitor change to 2,1 from 2
theta_min=theta(index);
err_squared=err^2;
end
