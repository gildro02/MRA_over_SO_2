% A function that take in all the relevant parameters to the 2D numerical
% experement (B, Q, N, sigma, the distribution rho, a, ex....) and outputs the
% resulting a and it's corresponding error (for both algorithms).
% ALSO OUTPUT THE MOMENT ERROR

function [error_squared_FM, coeff_FM, error_squared_spectral, coeff_spectral] = runNumericalExperiment(varargin)

% Create an input parser object
p = generateMRAInputParser();

% Parse the paramters
parse(p, varargin{:})
parameters = p.Results;

% Generate the moments, based on the type requested by the user:
if strcmp(parameters.moment_type, 'empirical') == 1
    [M_1, M_2] = generateEmpiricalMoments(parameters);
elseif strcmp(parameters.moment_type, 'analytical') == 1
    [M_1, M_2] = generateAnalyticalMoments(parameters);
else
    error("Illegal Input!");
end

% Run the algorithms
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