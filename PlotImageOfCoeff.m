% Input: a_symm_1B
% Output: figure with the plot of the corresponding picture to the
% parameters.
function fig = PlotImageOfCoeff(a_symm_1B, Phi_ns_mat, size_image, force_positive)

[image, region] = coeff2image(a_symm_1B, Phi_ns_mat, size_image);
if force_positive
    image = image - min(min(image, [], "all"), 0);
%     image(~region) = 0;
end

fig = figure;
imagesc(image, "AlphaData", double(region))
daspect([1 1 1])
set(gca, 'ydir', 'normal')
axis off
% colorbar
colormap gray
