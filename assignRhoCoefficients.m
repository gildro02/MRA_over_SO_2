% Input: Q, B, full coefficients of rho.
% Output: assigns all the relevent partial vectors of coefficients, and
% varrho.
function assignRhoCoefficients(Q, B, rho_hat_symm_2B)

% Warn of wrong use:
if inputname(1) == 'B'
    warning("Swap B and Q!");
end
%First, assign the frequency vectors:
assignFrequencies(Q, B);

% Define the vectors of interest:
rho_hat_symm_1B = rho_hat_symm_2B(-B <= k_symm_2B & k_symm_2B <= B);
rho_hat_half_2B = rho_hat_symm_2B(0 <= k_symm_2B);
rho_hat_half_1B = rho_hat_symm_2B(0 <= k_symm_2B & k_symm_2B <= B);
varrho = repelem(rho_hat_symm_1B, Q);

% Get the list of variables in the function's workspace to assign
variablesToAssign = {'rho_hat_symm_1B', 'rho_hat_half_2B', 'rho_hat_half_1B', 'varrho'};

% Assign those variables to the caller's workspace
for n = 1:length(variablesToAssign)
    assignin('caller', variablesToAssign{n}, eval(variablesToAssign{n}));
end

end