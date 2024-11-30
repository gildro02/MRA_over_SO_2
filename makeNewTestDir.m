% Input: -
% Output: makes a new directory in the archive of numerical test, adds 1 to
% test_number.
function [dir_path] = makeNewTestDir()

load('./Figures_Thesis/Comparison_2D_Archive/test_number')
test_number = test_number + 1;
save('./Figures_Thesis/Comparison_2D_Archive/test_number.mat','test_number')
test_end_string = "TestNumber_" + test_number;
mkdir('./Figures_Thesis/Comparison_2D_Archive', test_end_string);
dir_path = strcat('./Figures_Thesis/Comparison_2D_Archive/', test_end_string);
end