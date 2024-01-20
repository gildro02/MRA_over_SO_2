function [C] = circul(v)
%Create circulant matrix with first colomn v;
C=gallery("circul",v).';
end