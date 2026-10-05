outputPdf = fullfile(tempdir, 'Figure_2B_editable_source.pdf');
outputFig = ['/Users/xialeihuang/Desktop/Neuroscience_Project/' ...
    'Spectral_Analysis_Report/Figures/Final draft/Figure 2B.fig'];

pathwayNames = {'Leak', 'L6', 'Recurrent E', 'Inhibition', 'Net ODE'};

% Four angle-0 stimulus conditions from the original lower panel:
% contrasts 19, 42, 66, and 100, respectively.
conditionValues = [
    -1.0000, -1.0000, -1.0000, -1.0000
     0.9805,  1.0573,  0.9183,  0.9556
     0.8727,  0.8523,  0.6414,  0.6025
    -1.0686, -1.1292, -0.7621, -0.7337
    -0.2154, -0.2197, -0.2023, -0.1756
];

palette = [
    178  24  43
    227  74  51
    171 217 233
     67 147 195
     33 102 172
] / 255;

fig = figure('Color', 'w', 'Units', 'inches', ...
    'Position', [0.5 0.5 14.5 9.2], 'Renderer', 'painters');
ax = axes(fig, 'Position', [0.12 0.16 0.84 0.79]);
hold(ax, 'on');

violinHalfWidth = 0.24;
for pathwayIndex = 1:numel(pathwayNames)
    samples = conditionValues(pathwayIndex, :);
    sampleStd = std(samples, 0, 2);
    bandwidth = max(0.025, 0.55 * sampleStd);
    yGrid = linspace(min(samples) - 3 * bandwidth, ...
        max(samples) + 3 * bandwidth, 300);

    density = zeros(size(yGrid));
    for sampleIndex = 1:numel(samples)
        density = density + exp(-0.5 * ((yGrid - samples(sampleIndex)) / bandwidth).^2);
    end
    density = density / max(density);
    width = violinHalfWidth * density;

    patch(ax, [pathwayIndex - width, fliplr(pathwayIndex + width)], ...
        [yGrid, fliplr(yGrid)], palette(pathwayIndex, :), ...
        'FaceAlpha', 0.68, 'EdgeColor', palette(pathwayIndex, :), ...
        'LineWidth', 2.4);

end

yline(ax, 0, '-', 'Color', [0.35 0.35 0.35], ...
    'LineWidth', 1.6, 'HandleVisibility', 'off');

set(ax, 'FontName', 'Arial', 'FontSize', 26, 'LineWidth', 1.8, ...
    'Box', 'off', 'TickDir', 'out', 'XLim', [0.45 5.55], ...
    'XTick', 1:numel(pathwayNames), 'XTickLabel', pathwayNames, ...
    'YLim', [-1.50 1.30]);
ax.Toolbar.Visible = 'off';
ax.YGrid = 'on';
ax.XGrid = 'off';
ax.GridAlpha = 0.16;

ylabel(ax, 'Individual Contribution (1/\tau)', ...
    'Interpreter', 'tex', 'FontSize', 30, 'FontWeight', 'bold');

savefig(fig, outputFig);
exportgraphics(fig, outputPdf, 'ContentType', 'vector', ...
    'BackgroundColor', 'white');
