function [p] = generateMRAInputParser()

% Create an input parser object
p = inputParser;

% Define name-value pair parameters
% Syntax: addParameter(parserObject, 'Name', DefaultValue, ValidationFcn)

% Angular bandwidth
addParameter(p, 'B', [], @(x) length(x) == 1 & isnumeric(x) & (x > 0) & (mod(x, 1) == 0)) %check for positive whole number
% Radial nandwidth
addParameter(p, 'Q', [], @(x) length(x) == 1 & isnumeric(x) & (x > 0) & (mod(x, 1) == 0)) %check for positive whole number
% Number of model realiztions
addParameter(p, 'N', [], @(x) length(x) == 1 & isnumeric(x) & (x > 0) & (mod(x, 1) == 0)) %check for positive whole number
% Noise level
addParameter(p, 'sigma', [], @(x) length(x) == 1 & isnumeric(x) & (x >= 0)) %check for non-negative number
% Rotation distribution function
addParameter(p, 'rho_func', [], @(f) isa(f, 'function_handle')) %check for a function
% Fourier coefficients of rotation distribution function
dist_coeff_checker = @(coeff) isnumeric(coeff); %might add a more comperehensive condition.
addParameter(p, 'rho_hat_symm_2B', [], dist_coeff_checker)
% Vector of angles on which rho is a valid distribution (non-negative and
% integral = 1).
addParameter(p, 't', [], @(x) isreal(x)) % checks for a real vector
% The F-B coefficients of the original picture: (a vector of size
% (2*B+1)*Q)
addParameter(p, 'a_symm_1B', [], @(x) isnumeric(x));
% Binary option to use either empirical moments or analytical moments
moment_type_checker = @(x) (isstring(x) || ischar(x)) & (strcmp(x, 'empirical') || strcmp(x, 'analytical')); % checks for legal strings.
addParameter(p, 'moment_type', [], @(x) moment_type_checker(x));
% Resolution for the error determination:
% 1e-5 is about 1e-10 acccuracy and 0.75 seconds, 5e-4 is about 1e-5
% accuracy and 0.05 seconds (For FM with inf SNR).
addParameter(p, 'resolution_error', [], @(x) isreal(x) & x > 0 & length(x) == 1); %check for real positive number
% Options for only one algorithm
only_one_algorithm_checker = @(x) (isstring(x) || ischar(x)) & (strcmp(x, 'Frequency Marching') || strcmp(x, 'Spectral')); % checks for legal strings.
addParameter(p, 'only_one_algorithm', [], @(x) only_one_algorithm_checker(x)); %check for real positive number
% Forcing a_est_tilde to contain pure phases in the spectral algorithm
% Default - false
addParameter(p, 'force_pure_phases', false, @(x) islogical(x) & length(x) == 1) %check for true or false
end