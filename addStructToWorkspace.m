% Input: a struct of values, variable number of field names
% Output: adds all the values in the struct to the workspace by their
% names, only for those listed in the field names. If no field
% names are provided, all fields from the struct are added to the workspace.
function addStructToWorkspace(structData, varargin)

if isempty(varargin)
    % If no specific fields are mentioned, add all fields
    fieldNames = fieldnames(structData);
else
    % Otherwise, use the field names provided in varargin
    fieldNames = varargin;
end

% Loop through each field in fieldNames
for i = 1:length(fieldNames)
    fieldName = fieldNames{i};
    
    % Check if the field exists in the struct
    if isfield(structData, fieldName)
        % Get the field value and assign it to the workspace
        fieldValue = structData.(fieldName);
        assignin('caller', fieldName, fieldValue);
    else
        warning('Field "%s" does not exist in the struct.', fieldName);
    end
end
end
