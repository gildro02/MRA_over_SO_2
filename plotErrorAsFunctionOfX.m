% Input: the column vector of errors of the algorithms, and the vector of 
% corresponding parameters. The parameter vector can be of sigmas, N's, 
% M_2_diff or otherwise. The errors for the algorithms are relative to the 
% power of the image. Also, the number to average by is given and the 
% dist_from_circ is used for the plot.
function [fig] = plotErrorAsFunctionOfX(x_variable, x_name, error_squared_spectral_relative, error_squared_FM_relative, num_average, dist_from_circ)

% Check for legal input
if any(size(x_variable) ~= size(error_squared_spectral_relative)...
        | size(x_variable) ~= size(error_squared_FM_relative))
    error("Sizes do not match!");
end
if ~isstring(x_name) && ~ischar(x_name)
    error("Name is not a string/char vector!")
end
num_unique = length(x_variable) ./ num_average;
% If the plotting variable x is unsorted, sort it:
[x_sorted, I_sort] = sort(x_variable, "ascend");
error_squared_spectral_sorted = error_squared_spectral_relative(I_sort);
error_squared_FM_sorted = error_squared_FM_relative(I_sort);

% Take the mean of via blocks of size num_average:
x_average = mean(reshape(x_sorted, [num_average, num_unique]), 1).';
error_squared_spectral_mean = mean(reshape(error_squared_spectral_sorted,...
    [num_average, num_unique]), 1).';
error_squared_FM_mean = mean(reshape(error_squared_FM_sorted,...
    [num_average, num_unique]), 1).';

% plot Relative error:
fig = figure;
loglog(x_average, error_squared_spectral_mean);
hold on
loglog(x_average, error_squared_FM_mean);
xlabel(x_name, "fontsize", 15)
ylabel("Relative Squared Error", "fontsize", 15)

axis tight
legend("Spectral Algorithm", "Frequancy Marching algorithm");
title(sprintf("$||T-C||^2 = %.3f$", dist_from_circ)...
    , "fontsize", 15, "interpreter", "latex")


end