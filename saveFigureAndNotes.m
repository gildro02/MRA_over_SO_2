function saveFigureAndNotes(fig, fileName, folder, notes)
% saveFigureAndNotes Saves a figure and notes to a specified folder.
%   saveFigureAndNotes(fig, fileName, folder, notes) saves the figure
%   provided in `fig` as an image file with the name `fileName` in the
%   `folder` directory. Additionally, saves the string `notes` as a
%   text file in the same folder if provided.

% Validate inputs
if ~isgraphics(fig, 'figure')
    error('The first argument must be a valid figure handle.');
end
if ~ischar(fileName) && ~isstring(fileName)
    error('The second argument must be a valid file name as a string.');
end
if ~ischar(folder) && ~isstring(folder)
    error('The third argument must be a valid folder path as a string.');
end
if nargin < 4
    notes = ''; % Default to empty notes if not provided
end
if ~ischar(notes) && ~isstring(notes)
    error('The fourth argument must be a valid string for notes.');
end

% Ensure the folder exists
if ~exist(folder, 'dir')
    error('The specified folder does not exist.');
end

% Save the figure
figurePath = fullfile(folder, fileName);
saveas(fig, figurePath);

% Save the notes if provided
if ~isempty(notes)
    notesPath = fullfile(folder, 'notes.txt');
    fid = fopen(notesPath, 'w');
    if fid == -1
        error('Failed to create or open the notes file.');
    end
    fprintf(fid, '%s', notes);
    fclose(fid);
    fprintf('Notes saved as: %s\n', notesPath);
else
    fprintf('No notes provided; skipping notes file creation.\n');
end

fprintf('Figure saved as: %s\n', figurePath);
end
