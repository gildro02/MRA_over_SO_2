function [C] = BCCB(zeta)
% takes a Q x (2B+1) matrix zeta and returns the BCCB with zeta as its
% first column.
Q=size(zeta,1);
M=size(zeta,2);
B=(M-1)/2;
C_blocks_vec=cell(M,1);
if Q~=1
    for n=1:M
        C_blocks_vec{n}=circul(zeta(:,n)); %creates the QxQ circulant sub matrices
    end
    standard_circul=circul(1:M);
    C_blocks_mat=cell(M,M);
    for n=1:M
        for m=1:M
            C_blocks_mat{n,m}=C_blocks_vec{standard_circul(n,m)};
        end
    end
    C=cell2mat(C_blocks_mat);
else
    C=circul(zeta);
end