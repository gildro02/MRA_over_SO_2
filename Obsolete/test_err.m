%x,y are colomn vectors of length Q(2B+1)x1
M=2*B+1;
mat_x=mat(x,Q,B);
mat_y=mat(y,Q,B);
exp_multi_matrix=repmat(exp(1i.*(-B:B).*(2*pi/(2*B+1)).*permute(0:2*B,[1,3,2])),[Q,1,1]);
[err_2D,index]=min(vecnorm(reshape(exp_multi_matrix.*mat_x-mat_y,[Q*(2*B+1),1,2*B+1]),2,1),[],3);
l=index-1;
err_squared_2D=err_2D^2;

[err_1D,index]=min(vecnorm(x.*exp(1i.*(-B:B).'.*(2*pi)/(2*B+1).*(0:2*B))-y,2,1));
l=index-1;
err_squared_1D=err_1D^2;

reshape_1D=x.*exp(1i.*(-B:B).'.*(2*pi)/(2*B+1).*(0:2*B))-y;
reshape_2D=reshape(exp_multi_matrix.*mat_x-mat_y,[Q*(2*B+1),1,2*B+1]);

norm(squeeze(reshape_2D)-reshape_1D,"fro")^2