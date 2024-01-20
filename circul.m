% Input: v, a vector.
% Output: C, a circulant matrix with v as it's first column
function [C] = circul(v)
C=gallery("circul",v).';
end