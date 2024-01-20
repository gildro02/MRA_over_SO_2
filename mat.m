function [A] = mat(v,Q,B)
% returns the Qx2*B+1 matrix A of the Q(2*B+1)x1 vector v reorgenized to be Qx(2*B+1), going
% column by column (the inverse of vec);
A= reshape(v,[Q,2*B+1]);
end

