function [err_squared,l] = Circ_Error_Continuous_2D(x,y,Q,B)
%x,y are colomn vectors of length Q(2B+1)x1
M=2*B+1;
mat_x=mat(x,Q,B);
mat_y=mat(y,Q,B);
exp_multi_matrix=repmat(exp(1i.*(-B:B).*(2*pi/(2*B+1)).*permute(0:2*B,[1,3,2])),[Q,1,1]);
[err,index]=min(vecnorm(reshape(exp_multi_matrix.*mat_x-mat_y,[Q*(2*B+1),1,2*B+1]),2,1),[],3);
l=index-1;
err_squared=err^2;
end

