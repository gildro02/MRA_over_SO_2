% Input: Plotting parameters (variable for x axis, errors, number of
% reps,...), in name-value format.
% Output: the desired plot figure.
function [fig] = plotMRAError(varargin)

% Create an input parser object
p = inputParser;

% Define name-value pair parameters
% Syntax: addParameter(parserObject, 'Name', DefaultValue, ValidationFcn)

addParameter(p, 'variable', [], @(x) isreal(x) & isvector(x)) %check for real vector
addParameter(p, 'num_rep', [], @(x) length(x) == 1 & mod(x, 1) == 0 & (x > 0)) %check for natural number

mean_or_median_checker = @(x) (isstring(x) || ischar(x)) & (strcmp(x, 'mean') || strcmp(x, 'median'));
addParameter(p, 'mean_or_median', 'mean', @(x) mean_or_median_checker(x)) %checks for legal strings.
addParameter(p, 'error_spectral', [], @(x) isreal(x) & isvector(x)) %check for real vector
addParameter(p, 'error_FM', [], @(x) isreal(x) & isvector(x)) %check for real vector
addParameter(p, 'labels', [], @(x) iscell(x) & length(x) == 2) %check for 2-length cell array of xlabel and ylabel.
addParameter(p, 'title', "", @(x) isstring(x) || ischar(x)) %check for string/char vector.
addParameter(p, 'std_factor_or_prctile_region', [], @(x) isreal(x) & length(x) == 1 & x >= 0) %check for positive real number.
addParameter(p, 'is_x_math', false, @islogical) %check for logical value.
addParameter(p, 'y_asymptote_start_percent', [], @(x) isreal(x) & length(x) == 1 & x >= 0) %check for positive real number.
addParameter(p, 'y_asymptote_val', [], @(x) isreal(x) & length(x) == 1 & x >= 0) %check for positive real number.

% Parse the paramters
parse(p, varargin{:})
params = p.Results;

% Decide on mean or median of the errors:
if strcmp(params.mean_or_median, "mean")
    avg_func = @mean;
elseif strcmp(params.mean_or_median, "median")
    avg_func = @median;
else
    error("Specify mean or median!")
end

% Calculations before plot:
num_unique = length(params.variable) ./ params.num_rep;
% If the plotting variable x is unsorted, sort it:
[variable_sorted, I_sort] = sort(params.variable, "ascend");
error_spectral_sorted = params.error_spectral(I_sort);
error_FM_sorted = params.error_FM(I_sort);

% Take the mean/median of via blocks of size num_rep:
error_spectral_avg = avg_func(reshape(error_spectral_sorted,...
    [params.num_rep, num_unique]), 1).';
error_FM_avg = avg_func(reshape(error_FM_sorted,...
    [params.num_rep, num_unique]), 1).';

% If variable is not repetitions of a constant value, average each block:
variable_avg = mean(reshape(variable_sorted, [params.num_rep, num_unique]), 1).';

% plot the error:
fig = figure;
color_spectral = "blue";
color_FM = "red";

font_size = 15;
plot_spectral = loglog(variable_avg, error_spectral_avg, color_spectral, "LineWidth", 1.5);
hold on
plot_FM = loglog(variable_avg, error_FM_avg, color_FM, "LineWidth", 1.5);
if params.is_x_math 
    x_label_latex = strcat('$', params.labels{1}, '$');
else
    x_label_latex = params.labels{1};
end
y_label_latex = params.labels{2};
xlabel(x_label_latex, "fontsize", font_size, "Interpreter", "latex")
ylabel(y_label_latex, "fontsize", font_size, "Interpreter", "latex")

title_latex = params.title;
title(title_latex, "fontsize", font_size);

axis tight
grid on

% Plot semi-transperent region around the errors:
if ~isempty(params.std_factor_or_prctile_region)
    % Plot std region if mean is specified:
    if strcmp(params.mean_or_median, "mean")
        % Take the std of via blocks of size num_rep:
        std_spectral = std(reshape(error_spectral_sorted,...
            [params.num_rep, num_unique]), 0, 1).' .* params.std_factor_or_prctile_region;
        std_FM = std(reshape(error_FM_sorted,...
            [params.num_rep, num_unique]), 0, 1).' .* params.std_factor_or_prctile_region;
        
        std_spectral_upper = error_spectral_avg + std_spectral;
        std_spectral_lower = error_spectral_avg - std_spectral;
        
        std_FM_upper = error_FM_avg + std_FM;
        std_FM_lower = error_FM_avg - std_FM;
        
        % Plot the shaded area for the standard deviation
        fill([variable_avg; flip(variable_avg)], [std_spectral_upper; flip(std_spectral_lower)], ...
            color_spectral, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        
        fill([variable_avg; flip(variable_avg)], [std_FM_upper; flip(std_FM_lower)], ...
            color_FM, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    % Plot percentile region if median is specified:
    elseif  strcmp(params.mean_or_median, "median")
  
        prctile_spectral_upper = prctile(reshape(error_spectral_sorted,...
            [params.num_rep, num_unique]), 50 + params.std_factor_or_prctile_region, 1).';
        prctile_spectral_lower = prctile(reshape(error_spectral_sorted,...
            [params.num_rep, num_unique]), 50 - params.std_factor_or_prctile_region, 1).';
        
        prctile_FM_upper = prctile(reshape(error_FM_sorted,...
            [params.num_rep, num_unique]), 50 + params.std_factor_or_prctile_region, 1).';
        prctile_FM_lower = prctile(reshape(error_FM_sorted,...
            [params.num_rep, num_unique]), 50 - params.std_factor_or_prctile_region, 1).';
        
        % Plot the shaded area for the standard deviation
        fill([variable_avg; flip(variable_avg)], [prctile_spectral_upper; flip(prctile_spectral_lower)], ...
            color_spectral, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        
        fill([variable_avg; flip(variable_avg)], [prctile_FM_upper; flip(prctile_FM_lower)], ...
            color_FM, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    end
end

% Plot asymptote for y:
if ~isempty(params.y_asymptote_val)
    % Set Default start to the begining of the x axis.
    if isempty(params.y_asymptote_start_percent)
        params.y_asymptote_start_percent = 0;
    end
    x_limits = xlim;
    x_start = x_limits(1) .* (x_limits(2) ./ x_limits(1)) .^ ...
        (params.y_asymptote_start_percent / 100); % Start % away from the right edge
    x_end = x_limits(2);
    loglog([x_start, x_end], [params.y_asymptote_val, params.y_asymptote_val], 'k--', 'LineWidth', 1.5);
end
    
    
    
legend([plot_spectral, plot_FM], ...
    "Spectral Algorithm",...
    "Frequancy Marching",...
    "fontsize", font_size, "Interpreter", "latex");

