function make_CG_cluster_compare_6_8_10(outputPdf)
% Compare genuine CG-library E-population eigencluster envelopes.

artifactRoot = fileparts(mfilename('fullpath'));
if nargin < 1 || strlength(string(outputPdf)) == 0
    outputPdf = fullfile(artifactRoot, ...
        'CG_Library_Top_Eigenclusters_6_8_10.pdf');
end

loaded = load(fullfile(artifactRoot, 'cg_cluster_compare_6_8_10.mat'), ...
    'topEigenclusterEAll', 'clusterCounts');
maps = loaded.topEigenclusterEAll;
counts = loaded.clusterCounts;
assert(isequal(counts, [6 8 10]), 'Expected cluster counts [6 8 10].');

sharedLimits = [0, max(cellfun(@(x) max(x(:)), maps))];
fig = figure('Color', 'w', 'Units', 'inches', ...
    'Position', [0.5, 0.5, 18.5, 6.8]);
lefts = [0.055, 0.365, 0.675];
width = 0.245;
axesHandles = gobjects(1, 3);

for index = 1:3
    ax = axes(fig, 'Position', [lefts(index), 0.16, width, 0.72]);
    imagesc(ax, maps{index});
    axis(ax, 'image');
    set(ax, 'YDir', 'normal');
    colormap(ax, jet(256));
    caxis(ax, sharedLimits);
    xlim(ax, [0.5, 40.5]);
    ylim(ax, [0.5, 40.5]);
    xticks(ax, [1 10 20 30 40]);
    yticks(ax, [1 10 20 30 40]);
    xlabel(ax, 'Map column (pixel)', 'FontSize', 26);
    ylabel(ax, 'Map row (pixel)', 'FontSize', 26);
    title(ax, sprintf('Top %d Modes', counts(index)), ...
        'FontSize', 28, 'FontWeight', 'bold');
    set(ax, 'FontName', 'Arial', 'FontSize', 26, 'LineWidth', 1.25, ...
        'Box', 'on', 'TickDir', 'out', 'Layer', 'top');
    hold(ax, 'on');
    for boundary = [10.5 20.5 30.5]
        xline(ax, boundary, '-', 'Color', [0.45 0.45 0.45], ...
            'LineWidth', 0.8, 'HandleVisibility', 'off');
        yline(ax, boundary, '-', 'Color', [0.45 0.45 0.45], ...
            'LineWidth', 0.8, 'HandleVisibility', 'off');
    end
    hold(ax, 'off');
    axesHandles(index) = ax;
end

axisPosition = axesHandles(3).Position;
cb = colorbar(axesHandles(3), 'eastoutside');
axesHandles(3).Position = axisPosition;
cb.Position = [axisPosition(1) + axisPosition(3) + 0.008, ...
    axisPosition(2), 0.012, axisPosition(4)];
cb.FontName = 'Arial';
cb.FontSize = 26;
cb.LineWidth = 1.1;
cb.Title.String = 'sp/s';
cb.Title.FontName = 'Arial';
cb.Title.FontSize = 26;
cb.Ruler.Exponent = 0;

exportgraphics(fig, outputPdf, 'ContentType', 'image', 'Resolution', 600);
close(fig);
fprintf('Saved %s\n', outputPdf);
end
