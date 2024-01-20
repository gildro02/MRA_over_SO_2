function [Error,matched_y,opt_shift] = Circular_Error(x,y)
%x and y are the same length
%match second input to the first one
if size(x)~=size(y)
    error("non matching sizes")
end

L=length(x);
circ_error=zeros(1,L);
for shift=1:L
    circ_error(shift)=norm(x-circshift(y,shift),"fro").^2;
end
[Error,opt_shift]=min(circ_error);
matched_y=circshift(y,opt_shift);
