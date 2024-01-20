function [C] = circul_semi(v,L)
%v=first colomn of semi circulant matrix
%L=number of colomns. L must divide length(v), and r=length(v)/L is the
%jump in indecies between each colomn.
if size(v,1)==1
    v=v.';
end
r=length(v)/L;
if r~=floor(r)
    error("length(v)/L is not an integer");
end
C=zeros(r*L,L);
for n=0:L-1
    C(:,n+1)=circshift(v,r*n);
end

end

