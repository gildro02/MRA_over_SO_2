% Input: 1) coeff, a column vector of the F-B coefficients of the image
% 2) Phi_ns_mat, a matrix with columns that are the basis vectors, ordered
% like coeff.
% 3) B, Q.
% Output: a matrix representing the image.
function [image] = coeff2image(coeff, Phi_ns_mat, size_image)
%% reshape coeff to be a column vector:
coeff=reshape(coeff,[length(coeff) 1]);
%%
image_vec=Phi_ns_mat*coeff;
size_half=floor(size_image/2);
[x, y]=meshgrid(-size_half:size_half, -size_half:size_half);
r=sqrt(x.^2+y.^2);
image=zeros(size_image,size_image);
image(r<=size_half) = image_vec;
if norm(imag(image),"fro")>1e-10
    warning("Picture is not Real!");
    %error("Picture is not Real!");
end
image=real(image); % ensure realness
end