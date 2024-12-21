% Input: the parameters of the algorithm, optional - specify only 'Spectral' or
% 'Frequency Marching' Algorithms.
% Output: The output coefficients of the algorithms and the corresponding
% errors, or [] for the case where only the other algorithm was selected.
function [error_squared_FM, coeff_FM, error_squared_spectral, coeff_spectral] = methodOfMomentsAlgorithmsHandler(M_1, M_2, parameters)

if strcmp(parameters.only_one_algorithm, 'Spectral') == 1 % If specified to only run the spectral algorithm
    [coeff_spectral, error_squared_spectral] = spectralAlgorithm(M_1, M_2, parameters);
    coeff_FM = [];
    error_squared_FM = [];
elseif strcmp(parameters.only_one_algorithm, 'Frequency Marching') ==  1 % If specified to only run the FM algorithm
    [coeff_FM, error_squared_FM] = frequencyMarchingAlgorithm(M_1, M_2, parameters);
    coeff_spectral = [];
    error_squared_spectral = [];
elseif isempty(parameters.only_one_algorithm) == 1 % If nothing specified (default) - run both
    [coeff_spectral, error_squared_spectral] = spectralAlgorithm(M_1, M_2, parameters);
    [coeff_FM, error_squared_FM] = frequencyMarchingAlgorithm(M_1, M_2, parameters);
else % Error in specification
    error("Illegal Input!");
end


end