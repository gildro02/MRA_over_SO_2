% Input: zeta, a Q x (2B+1) matrix zeta
% Output: BTTB with vec(zeta) as its first column
function [T] = BTTB(zeta)
Q=size(zeta,1);
M=size(zeta,2);
B=(M-1)/2;
T_blocks_vec=cell(M,1);
for col=1:M
    T_blocks_vec{col}=toeplitz(zeta(:,col),zeta(:,col)); %creates the QxQ toeplitz sub matrices
end
standard_toeplitz=toeplitz(1:M);
T_blocks_mat=cell(M,M);
for row=1:M
    for col=1:M
        if col<row
            T_blocks_mat{row,col}=T_blocks_vec{standard_toeplitz(row,col)};
        else
            T_blocks_mat{row,col}=T_blocks_vec{standard_toeplitz(row,col)}';
        end
    end
end
T=cell2mat(T_blocks_mat);

