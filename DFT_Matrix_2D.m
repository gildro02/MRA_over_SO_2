function [W] = DFT_Matrix_2D(Q, B)
M = 2 * B + 1;
W_M = (1 / sqrt(M)) * dftmtx(M);
W_Q = (1 / sqrt(Q)) * dftmtx(Q);
W = kron(W_M, W_Q);
end