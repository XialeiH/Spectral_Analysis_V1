function make_Final_Figure_1B
% Export the left panel of Figure 1C without a title.

artifactRoot = fileparts(mfilename('fullpath'));
sourceFig = fullfile(artifactRoot, ...
    '1.6_Paper3_CG_library_vs_h96_DNN_pixelwise_firing_rates_10pct.fig');
outputDir = ['/Users/xialeihuang/Desktop/Neuroscience_Project/' ...
    'Spectral_Analysis_Report/Figures/Final draft'];
outputPdf = fullfile(outputDir, 'Figure 1B.pdf');
outputFig = fullfile(outputDir, 'Figure 1B.fig');

neuroColors = [
    178  24  43
    227  74  51
    254 227 145
    171 217 233
     67 147 195
     33 102 172
] / 255;

sourceFigure = openfig(sourceFig, 'invisible');
sourceAxes = findall(sourceFigure, 'Type', 'axes');
[~, order] = sort(arrayfun(@(sourceAx) sourceAx.Position(1), sourceAxes));
sourceAxes = sourceAxes(order);
assert(numel(sourceAxes) == 2, ...
    'Expected exactly two axes in the Figure 1C source.');

fig = figure('Visible', 'off', 'Units', 'inches', ...
    'Position', [0.5 0.5 9.5 8.5], 'Color', 'w', ...
    'Renderer', 'painters');
ax = axes(fig, 'Position', [0.165 0.145 0.800 0.815]);
copyobj(allchild(sourceAxes(1)), ax);
close(sourceFigure);
ax.Toolbar.Visible = 'off';
set(ax, 'FontName', 'Arial', 'FontSize', 26, 'LineWidth', 1.6, ...
    'Box', 'off', 'TickDir', 'out', 'Layer', 'top');

scatterObject = findobj(ax, 'Type', 'scatter');
assert(isscalar(scatterObject), 'Expected one scatter object.');
set(scatterObject, 'MarkerEdgeColor', neuroColors(1,:), ...
    'MarkerFaceColor', 'none');

title(ax, '');
xlabel(ax, 'CG Firing Rate (sp/s)', 'FontName', 'Arial', ...
    'FontSize', 26, 'FontWeight', 'normal');
ylabel(ax, 'DNN Surrogate Firing Rate (sp/s)', 'FontName', 'Arial', ...
    'FontSize', 26, 'FontWeight', 'normal');
xlim(ax, [3 40]);
ylim(ax, [3 40]);
xticks(ax, [3 10 20 30 40]);
yticks(ax, [3 10 20 30 40]);

savefig(fig, outputFig);
exportgraphics(fig, outputPdf, 'ContentType', 'vector', ...
    'BackgroundColor', 'white');
close(fig);
fprintf('Saved %s and %s.\n', outputPdf, outputFig);
end
