function make_Final_Figure_3D()
% Figure 3D.3: count distributions of log10 singular values.

scriptDir = fileparts(mfilename('fullpath'));
dataFile = fullfile(scriptDir, 'figure_3d3_log_sigma.tsv');
outputPdf = fullfile(tempdir, 'Figure_3D_editable_source.pdf');
outputFig = ['/Users/xialeihuang/Desktop/Neuroscience_Project/' ...
    'Spectral_Analysis_Report/Figures/Final draft/Figure 3D.fig'];

T = readtable(dataFile, 'FileType', 'text', 'Delimiter', '\t', ...
    'VariableNamingRule', 'preserve');
panelNames = {'Baseline Gaussian Connectivity', ...
    'Row-Sum-Matched Flattened Connectivity'};

neuroColors = [
    178  24  43
    227  74  51
    254 227 145
    171 217 233
     67 147 195
     33 102 172
] / 255;
lowColor = neuroColors(1, :);
highColor = neuroColors(6, :);

binWidth = 0.05;
xLimits = [-3.85 1.45];
binEdges = xLimits(1):binWidth:xLimits(2);
binCenters = binEdges(1:end-1) + binWidth / 2;
threshold = log10(0.05);

allCounts = cell(1, 2);
for panelIdx = 1:2
    mask = strcmp(T.panel, panelNames{panelIdx});
    allCounts{panelIdx} = histcounts(T.log10_sigma(mask), binEdges);
    assert(sum(allCounts{panelIdx}) == 4800, ...
        'Each panel must contain 4,800 singular values.');
end
yMaximum = 1.08 * max(cellfun(@max, allCounts));

fig = figure('Color', 'w', 'Units', 'inches', ...
    'Position', [0.5 0.5 18.0 7.2]);
tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

for panelIdx = 1:2
    ax = nexttile(tl, panelIdx);
    hold(ax, 'on');
    counts = allCounts{panelIdx};
    lowMask = binCenters <= threshold;
    highMask = ~lowMask;

    bar(ax, binCenters(lowMask), counts(lowMask), 1.0, ...
        'FaceColor', lowColor, 'EdgeColor', 'none');
    bar(ax, binCenters(highMask), counts(highMask), 1.0, ...
        'FaceColor', highColor, 'EdgeColor', 'none');
    xline(ax, threshold, '--', 'Color', [0.30 0.30 0.30], ...
        'LineWidth', 2.2, 'HandleVisibility', 'off');

    xlim(ax, xLimits);
    ylim(ax, [0 yMaximum]);
    xticks(ax, -3:1:1);
    xlabel(ax, '$\log_{10}(\sigma_j)$', 'Interpreter', 'latex');
    ylabel(ax, 'Count');
    title(ax, panelNames{panelIdx}, 'FontWeight', 'bold');
    set(ax, 'FontName', 'Arial', 'FontSize', 26, 'LineWidth', 1.5, ...
        'TickDir', 'out', 'Box', 'off', 'Layer', 'top');
    grid(ax, 'off');

    if panelIdx == 2
        legend(ax, {'$\sigma_j \leq 0.05$', '$\sigma_j > 0.05$'}, ...
            'Interpreter', 'latex', 'Location', 'northwest', ...
            'FontSize', 24, 'Box', 'off');
    end
end

savefig(fig, outputFig);
exportgraphics(fig, outputPdf, 'ContentType', 'vector');
fprintf('Saved %s\n', outputPdf);
end
