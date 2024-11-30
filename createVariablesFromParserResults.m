% Input: input parser struct p.
% Output: Assigns the values in p to the workspace of the caller. Throws an
% error if a parameter field is empty.

function createVariablesFromParserResults(p)
% Extract the field names from the Results struct
fields = fieldnames(p.Results);

% Dynamically create variables in the caller's workspace
for n = 1:numel(fields)
    if isempty(p.Results.(fields{n}))
        error(fields{n} + " is empty!")
    else
        assignin('caller', fields{n}, p.Results.(fields{n}));
    end
end
end