function saveFigureToAllFormats(fig, fileName, folder, notes)
% saveFigureToAllFormats Saves a figure and notes to a specified folder.
%   saveFigureToAllFormats(fig, fileName, folder, notes) saves the figure
%   provided in `fig` as an image file with the name `fileName` in the
%   `folder` directory. It is saved as a .fig, .eps and .pdf file.
%   Additionally, saves the string `notes` as a text file in the same 
%   folder if provided.

if nargin < 4
    saveFigureAndNotes(fig, strcat(fileName, '.fig'), folder)
else
    saveFigureAndNotes(fig, strcat(fileName, '.fig'), folder, notes);
end
saveFigureToEPS(fig, strcat(fileName, '.eps'), folder);
saveFigureToPDF(fig, strcat(fileName, '.pdf'), folder);
end