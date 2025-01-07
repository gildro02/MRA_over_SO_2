function saveFigureToPDF(figHandle, fileName, destFolder)
    % saveFigureToPDF - Save a MATLAB figure to a specified folder in PDF format
    %
    % Syntax:
    %   saveFigureToPDF(figHandle, fileName, destFolder)
    %
    % Inputs:
    %   figHandle  - (Optional) Handle of the figure to save. If empty or not provided, uses the current figure (gcf).
    %   fileName   - Name of the PDF file (with or without .pdf extension).
    %   destFolder - Relative path to the destination folder.
    %
    % Example:
    %   saveFigureToPDF(gcf, 'my_figure', 'output')
    %

    % Check if the figure handle is provided; otherwise, use the current figure
    if isempty(figHandle)
        figHandle = gcf;
    end

    % Validate the destination folder and ensure it exists
    if isempty(destFolder)
        error('Destination folder must be specified.');
    end

    % Ensure the file name has the .pdf extension
    if isempty(fileName)
        error('File name must be specified.');
    end
    [~, ~, ext] = fileparts(fileName);
    if isempty(ext)
        fileName = [fileName, '.pdf'];
    elseif ~strcmpi(ext, '.pdf')
        error('File name must have a .pdf extension.');
    end

    % Construct the full file path
    fullFilePath = fullfile(destFolder, fileName);

    % Save the figure in PDF format
    print(figHandle, fullFilePath, '-dpdf', '-r1000'); % '-r1000' for high resolution

    % Notify the user
    fprintf('Figure saved as PDF to: %s\n', fullFilePath);
end
