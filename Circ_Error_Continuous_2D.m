% Input: x,y are column vectors of the same length Q(2B+1)
% Output: 1) The circular squared error between them, assuming that they represent the
% Fourier-Bessel coefficients of a 2D function, defined on the unit disk. This rotation
% is restricted to quantization by steps of size 2pi/(2B+1).
% 2) l, the index of rotation of the first vector where the mininum
% distance is met (and the corresponding angle of rotation is
% phi = l*2pi/(2B+1)).
function [err_squared,l] = Circ_Error_Continuous_2D(x,y,Q,B)

M=2*B+1;
mat_x=mat(x,Q,B);
mat_y=mat(y,Q,B);
exp_multi_matrix=repmat(exp(1i.*(-B:B).*(2*pi/(2*B+1)).*permute(0:2*B,[1,3,2])),[Q,1,1]);
[err,index]=min(vecnorm(reshape(exp_multi_matrix.*mat_x-mat_y,[Q*(2*B+1),1,2*B+1]),2,1),[],3);
l=index-1;
err_squared=err^2;
end

