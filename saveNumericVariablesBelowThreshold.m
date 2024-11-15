function saveNumericVariablesBelowThreshold(folderPath, fileName, sizeThresholdKB)
    % Set the size threshold in bytes (1 KB = 1024 bytes)
    sizeThresholdBytes = sizeThresholdKB * 1024;

    % Get all variables in the base workspace
    variablesInfo = evalin('base', 'whos');
    
    % Construct the full file path
    filePath = fullfile(folderPath, strcat(fileName, '.mat'));

    % Initialize an empty structure to store variables that meet the criteria
    dataStruct = struct();

    % Loop through all variables and filter numeric ones below the threshold
    for i = 1:length(variablesInfo)
        varName = variablesInfo(i).name;
        varSize = variablesInfo(i).bytes;
        varClass = variablesInfo(i).class;
        
        % Check if the variable is numeric and below the size threshold
        if isnumeric(evalin('base', varName)) && varSize < sizeThresholdBytes
            dataStruct.(varName) = evalin('base', varName);
        end
    end

    % Save the filtered variables in the specified .mat file
    save(filePath, '-struct', 'dataStruct');
    
    fprintf('Numeric variables below %d KB saved to %s\n', sizeThresholdKB, filePath);
end
