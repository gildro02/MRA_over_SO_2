% A function that take in all the relevant parameters to the 2D numerical
% experement (B, Q, N, sigma, the distribution rho, a, ex....) and outputs the
% resulting a and it's corresponding error (for both algorithms).
% Also, if "empirical" moments are specified, outputs the frovenious error
% (squared) between the analytical and empirical moments.

function [error_squared_FM, coeff_FM, error_squared_spectral, coeff_spectral, M_1_diff, M_2_diff] = runNumericalExperiment(varargin)

% Create an input parser object
p = generateMRAInputParser();

% Parse the paramters
parse(p, varargin{:})
parameters = p.Results;

% Generate the moments, based on the type requested by the user:
[M_1, M_1_diff, M_2, M_2_diff] = momentHandler(parameters);

% Run the algorithms
[error_squared_FM, coeff_FM, error_squared_spectral, coeff_spectral] = ...
    methodOfMomentsAlgorithmsHandler(M_1, M_2, parameters);

end