% Input: x,y are column vectors of the same length 2B+1
% Output: 1) The circular squared error between them, assuming that they represent the
% Fourier coefficients of a 1D function, defined on [0,2pi]. This rotation
% is restricted to quantization by steps of size 2pi/(2B+1).
% 2) l, the index of rotation of the first vector where the mininum
% distance is met (and the corresponding angle of rotation is
% phi = l*2pi/(2B+1)).
function [err_squared,l] = circ_error_continuous(x,y,B)
[err,index]=min(vecnorm(x.*exp(1i.*(-B:B).'.*(2*pi)/(2*B+1).*(0:2*B))-y,2));
l=index-1;
err_squared=err^2;
end

