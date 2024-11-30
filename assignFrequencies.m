% Input: Q,B
% Output: Adds all the relevant frequency vectors and matricies to the
% context (workspace) of the caller (exept Q and B which already exist).
function assignFrequencies(Q,B)
% Warn of wrong use:
if inputname(1) == 'B'
    warning("Swap B and Q!");
end

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

% Get the list of variables in the function's workspace
allVariables = who;

% Get the names of function input arguments
varNames = arrayfun(@inputname, 1:nargin, 'UniformOutput', false);

% Filter out parameters from the list of variables
% This is as to not add Q and B again.
variablesToAssign = setdiff(allVariables, varNames);

% Assign remaining variables to the caller's workspace
for n = 1:length(variablesToAssign)
    assignin('caller', variablesToAssign{n}, eval(variablesToAssign{n}));
end
end