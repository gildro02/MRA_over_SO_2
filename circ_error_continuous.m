function [err_squared,l] = circ_error_continuous(x,y,B)
%x,y are colomn vectors.
[err,index]=min(vecnorm(x.*exp(1i.*(-B:B).'.*(2*pi)/(2*B+1).*(0:2*B))-y,2));
l=index-1;
err_squared=err^2;
end

