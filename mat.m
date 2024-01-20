% Input: v, a vector of length Q(2B+1).
% Output: A, a matrix of size Qx2*B+1 of the vector v, reorgenized to be Qx(2*B+1), going
% column by column (the inverse of vec);
function [A] = mat(v,Q,B)
A= reshape(v,[Q,2*B+1]);
end

