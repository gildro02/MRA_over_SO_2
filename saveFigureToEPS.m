function saveFigureToEPS(figHandle, fileName, destFolder)
    % saveFigureToEPS - Save a MATLAB figure to a specified folder in EPS format
    %
    % Syntax:
    %   saveFigureToEPS(figHandle, destFolder, fileName)
    %
    % Inputs:
    %   figHandle  - (Optional) Handle of the figure to save. If empty or not provided, uses the current figure (gcf).
    %   fileName   - Name of the EPS file (with or without .eps extension).
    %   destFolder - Relative path to the destination folder.
    %  
    %
    % Example:
    %   saveFigureToEPS(gcf, 'my_figure', 'output')
    %

    % Check if the figure handle is provided; otherwise, use the current figure
    if isempty(figHandle)
        figHandle = gcf;
    end

    % Validate the destination folder and ensure it exists
    if isempty(destFolder)
        error('Destination folder must be specified.');
    end

    % Ensure the file name has the .eps extension
    if isempty(fileName)
        error('File name must be specified.');
    end
    [~, ~, ext] = fileparts(fileName);
    if isempty(ext)
        fileName = [fileName, '.eps'];
    elseif ~strcmpi(ext, '.eps')
        error('File name must have an .eps extension.');
    end

    % Construct the full file path
    fullFilePath = fullfile(destFolder, fileName);

    % Save the figure in EPS format
    print(figHandle, fullFilePath, '-depsc', '-r1000'); % '-r1000' for high resolution

    % Notify the user
    fprintf('Figure saved as EPS to: %s\n', fullFilePath);
end
