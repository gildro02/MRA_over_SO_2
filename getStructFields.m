function varargout = getStructFields(s, varargin)
    % getStructFields: Retrieve values of specified struct fields.
    %
    % Usage:
    %   [v1, v2, ...] = getStructFields(s, 'field1', 'field2', ...)
    %
    % Input:
    %   s: Struct containing the fields.
    %   varargin: Variable number of character vectors representing field names.
    %
    % Output:
    %   varargout: Values of the requested fields in the specified order.
    %              Returns [] for non-existing fields.

    % Validate input
    if ~isstruct(s)
        error('First input must be a struct.');
    end
    
    % Preallocate output
    nFieldsRequested = length(varargin);
    varargout = cell(1, nFieldsRequested);
    
    % Retrieve values for each requested field
    for k = 1:nFieldsRequested
        fieldName = varargin{k};
        if isfield(s, fieldName)
            varargout{k} = s.(fieldName); % Assign field value
        else
            varargout{k} = []; % Assign empty for non-existent field
        end
    end
end
