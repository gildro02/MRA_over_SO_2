function [file_name] = saveFrequencies(Q,B)
% Saves relevant frequency vector
k_symm_1B=-B:B;
k_half_1B=0:B;
k_symm_2B=-2*B:2*B;
k_half_2B=0:2*B;

k_symm_1B_Q_mat=repmat(k_symm_1B,[Q,1]);
k_half_1B_Q_mat=repmat(k_half_1B,[Q,1]);
k_symm_2B_Q_mat=repmat(k_symm_2B,[Q,1]);
k_half_2B_Q_mat=repmat(k_half_2B,[Q,1]);

k_symm_1B_Q_vec=vec(k_symm_1B_Q_mat);
k_half_1B_Q_vec=vec(k_half_1B_Q_mat);
k_symm_2B_Q_vec=vec(k_symm_2B_Q_mat);
k_half_2B_Q_vec=vec(k_half_2B_Q_mat);

q_1B=0:Q-1;
q_full_1B_mat=repmat(q_1B.',[1,2*B+1]);
q_full_1B_vec=vec(q_full_1B_mat);

file_name='freq.mat';
save(file_name,'k_symm_1B','k_half_1B','k_symm_2B','k_half_2B',...
    'k_symm_1B_Q_mat', 'k_half_1B_Q_mat', 'k_symm_2B_Q_mat', 'k_half_2B_Q_mat',...
    'k_symm_1B_Q_vec', 'k_half_1B_Q_vec', 'k_symm_2B_Q_vec', 'k_half_2B_Q_vec',...
    'q_1B','q_full_1B_mat', 'q_full_1B_vec','-mat');