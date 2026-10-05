function make_Final_Figure_1C
% Combine the two pixel schematics and their E tuning curves.

artifactRoot = fileparts(mfilename('fullpath'));
sourceMat = fullfile(artifactRoot, ...
    '1.7_H96_E_orientation_tuning_curves_contrasts_19_42_66_100_data.mat');
outputDir = ['/Users/xialeihuang/Desktop/Neuroscience_Project/' ...
    'Spectral_Analysis_Report/Figures/Final draft'];
outputPdf = fullfile(outputDir, 'Figure 1C.pdf');
outputFig = fullfile(outputDir, 'Figure 1C.fig');

loaded = load(sourceMat, 'E_curves', 'provenance');
eCurves = loaded.E_curves;
sourceAngles = loaded.provenance.extendedAngles(:);
contrasts = loaded.provenance.contrasts(:).';
eCurves = eCurves(sourceAngles < 180, :, :);

neuroColors = [
    178  24  43
    227  74  51
    254 227 145
    171 217 233
     67 147 195
     33 102 172
] / 255;
curveColors = neuroColors([4 5 2 1], :);
xWindows = {-90:3.75:90, -67.5:3.75:112.5};
xTicks = {[-90 -45 0 45 90], [-67.5 -22.5 22.5 67.5 112.5]};

fig = figure('Visible', 'off', 'Color', 'w', 'Units', 'inches', ...
    'Position', [0.4 0.4 19 12], 'Renderer', 'painters');

schematicPositions = [0.070 0.720 0.160 0.253; ...
                      0.550 0.720 0.160 0.253];
curvePositions = [0.070 0.085 0.400 0.500; ...
                  0.550 0.085 0.400 0.500];
pixelRows = [5 1];

for panelIndex = 1:2
    schematicAxes = axes(fig, 'Position', ...
        schematicPositions(panelIndex, :));
    draw_pixel_panel(schematicAxes, pixelRows(panelIndex), 1, panelIndex);

    curveAxes = axes(fig, 'Position', curvePositions(panelIndex, :));
    hold(curveAxes, 'on');
    plotAngles = xWindows{panelIndex};
    wrappedAngles = mod(plotAngles, 180);
    sampleIndex = round(wrappedAngles / 3.75) + 1;
    curveLines = gobjects(1, numel(contrasts));
    for contrastIndex = 1:numel(contrasts)
        curve = eCurves(sampleIndex, contrastIndex, panelIndex);
        curveLines(contrastIndex) = plot(curveAxes, plotAngles, curve, '-', ...
            'Color', curveColors(contrastIndex, :), ...
            'LineWidth', 3.5, 'Marker', 'none');
    end
    xlim(curveAxes, [plotAngles(1) plotAngles(end)]);
    xticks(curveAxes, xTicks{panelIndex});
    ylim(curveAxes, [0 32]);
    yticks(curveAxes, 0:10:30);
    xlabel(curveAxes, 'Orientation (deg)', 'FontName', 'Arial', ...
        'FontSize', 26, 'FontWeight', 'normal');
    ylabel(curveAxes, 'E Firing Rate (sp/s)', 'FontName', 'Arial', ...
        'FontSize', 26, 'FontWeight', 'normal');
    grid(curveAxes, 'on');
    set(curveAxes, 'FontName', 'Arial', 'FontSize', 26, ...
        'LineWidth', 1.6, 'Box', 'off', 'TickDir', 'out', ...
        'Layer', 'top', 'GridAlpha', 0.18, 'MinorGridAlpha', 0.10);
    curveAxes.Toolbar.Visible = 'off';

    if panelIndex == 2
        legendLabels = arrayfun(@(contrast) sprintf('%d', contrast), ...
            contrasts, 'UniformOutput', false);
        legendHandle = legend(curveAxes, curveLines, legendLabels, ...
            'Location', 'northwest', 'NumColumns', 2, ...
            'FontName', 'Arial', 'FontSize', 24, 'Box', 'on');
        legendHandle.Title.String = 'Contrast';
        legendHandle.Title.FontSize = 24;
        legendHandle.Color = 'w';
        legendHandle.EdgeColor = [0.25 0.25 0.25];
        legendHandle.LineWidth = 1.0;
    end
end

savefig(fig, outputFig);
exportgraphics(fig, outputPdf, 'ContentType', 'vector', ...
    'BackgroundColor', 'white');
close(fig);
fprintf('Saved %s and %s.\n', outputPdf, outputFig);
end

function draw_pixel_panel(ax, pixelRow, pixelColumn, pixelLabel)
gridColor = [0.30 1.00 0.62];
pixelColor = [0.05 1.00 0.00];
diagonalColor = [0.48 0.48 1.00];
barColor = [1.00 0.08 0.08];

hold(ax, 'on');
axis(ax, 'equal');
axis(ax, [0 14 0 14]);
axis(ax, 'off');
ax.Clipping = 'off';
ax.Toolbar.Visible = 'off';

x0 = 2 + pixelColumn - 1;
y0 = 12 - pixelRow;
for coordinate = 0:14
    plot(ax, [coordinate coordinate], [0 14], '-', ...
        'Color', gridColor, 'LineWidth', 2.2);
    plot(ax, [0 14], [coordinate coordinate], '-', ...
        'Color', gridColor, 'LineWidth', 2.2);
end

draw_boundaries(ax);
guides = {
    [0 2], [0 2]; [0 2], [2 0]; ...
    [12 14], [0 2]; [12 14], [2 0]; ...
    [0 2], [12 14]; [0 2], [14 12]; ...
    [12 14], [12 14]; [12 14], [14 12]; ...
    [2 12], [2 12]; [2 12], [12 2]};
for guideIndex = 1:size(guides, 1)
    plot(ax, guides{guideIndex, 1}, guides{guideIndex, 2}, '--', ...
        'Color', diagonalColor, 'LineWidth', 2.4);
end

scatter(ax, [2 12 2 12], [2 2 12 12], 85, ...
    'MarkerFaceColor', [0 0 0.55], 'MarkerEdgeColor', 'none');
plot(ax, [6.1 7.9], [10.9 9.1], '-', ...
    'Color', barColor, 'LineWidth', 7);
plot(ax, [3.5 3.5], [5.7 8.3], '-', ...
    'Color', barColor, 'LineWidth', 7);
plot(ax, [8.2 10.8], [7.0 7.0], '-', ...
    'Color', barColor, 'LineWidth', 7);
plot(ax, [6.1 7.9], [3.1 4.9], '-', ...
    'Color', barColor, 'LineWidth', 7);

rectangle(ax, 'Position', [x0 y0 1 1], ...
    'FaceColor', pixelColor, 'EdgeColor', pixelColor, 'LineWidth', 3.0);
draw_boundaries(ax);
scatter(ax, [2 12 2 12], [2 2 12 12], 85, ...
    'MarkerFaceColor', [0 0 0.55], 'MarkerEdgeColor', 'none');
text(ax, x0 + 0.5, y0 + 0.5, sprintf('%d', pixelLabel), ...
    'FontName', 'Arial', 'FontSize', 28, 'FontWeight', 'bold', ...
    'Color', 'k', 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', 'Clipping', 'off');
end

function draw_boundaries(ax)
plot(ax, [2 2], [0 14], 'k-', 'LineWidth', 2.5);
plot(ax, [12 12], [0 14], 'k-', 'LineWidth', 2.5);
plot(ax, [0 14], [2 2], 'k-', 'LineWidth', 2.5);
plot(ax, [0 14], [12 12], 'k-', 'LineWidth', 2.5);
end
