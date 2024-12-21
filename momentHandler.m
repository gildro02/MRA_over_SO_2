% Input: the parameters of the algorithm, option 'empirical' or
% 'analytical' for moment type.
% Output: the correct moments, and, for the empirical case, the distance of
% the moment from the true moment (which is zero for the analytical case).
function [M_1, M_1_diff, M_2, M_2_diff] = momentHandler(parameters)

% Generate the analytical moments:
[M_1_analytical, M_2_analytical] = generateAnalyticalMoments(parameters);

% Empirical case:
if strcmp(parameters.moment_type, 'empirical') == 1
    [M_1, M_2] = generateEmpiricalMoments(parameters);
    M_1_diff = norm(M_1 - M_1_analytical, "fro") .^ 2;
    M_2_diff = norm(M_2 - M_2_analytical, "fro") .^ 2;
% Analytical case:
elseif strcmp(parameters.moment_type, 'analytical') == 1
    M_1 = M_1_analytical;
    M_2 = M_2_analytical;
    M_1_diff = 0;
    M_2_diff = 0;
else
    error("Illegal Input!");
end
end