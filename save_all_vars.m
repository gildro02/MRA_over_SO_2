% Your testing code here to create variables

% Get all variables in the workspace
vars = whos;

% Define the threshold for exclusion (1 million doubles)
size_threshold = 1e6;

% Create a structure to store variables
data = struct();

% Loop through all variables and store them in the structure, excluding large ones
for i = 1:numel(vars)
    var_name = vars(i).name;
    var_size = vars(i).bytes / 8; % Convert bytes to doubles
    
    % Check if the variable size is below the threshold
    if var_size <= size_threshold
        % Add "1D" suffix to the variable name
        new_var_name = [var_name '_2D'];
        
        % Store the variable in the structure with the new name
        data.(new_var_name) = eval(var_name);
    else
        disp(['Variable ' var_name ' exceeds the size threshold and will be excluded.']);
    end
end

% Save the structure to a .mat file
save('all_variables_2D.mat', '-struct', 'data');
