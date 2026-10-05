function make_Final_Figure_4B
% Plot normalized I right and left singular modes for J_full and J_I.

artifactRoot = fileparts(mfilename('fullpath'));
dataFile = fullfile(artifactRoot, 'Figure_4B_2_singular_I_maps.mat');
outputPdf = fullfile(tempdir, 'Figure_4B_editable_source.pdf');
outputFig = ['/Users/xialeihuang/Desktop/Neuroscience_Project/' ...
    'Spectral_Analysis_Report/Figures/Final draft/Figure 4B.fig'];
if ~isfile(dataFile)
    prepare_Figure_4B_2_I_maps;
end
loaded = load(dataFile, 'singularMaps');

rowFields = {'J_full', 'J_I'};
columnFields = {'RightSingular', 'LeftSingular'};
columnTitles = {'Right Singular Mode', 'Left Singular Mode'};
jacobianLabels = {'J_{full}', 'J_I'};
columnLimits = zeros(2, 2);
for columnIndex = 1:2
    maximum = 0;
    for rowIndex = 1:2
        values = loaded.singularMaps.(rowFields{rowIndex}). ...
            (columnFields{columnIndex});
        maximum = max(maximum, max(abs(values), [], 'all'));
    end
    columnLimits(columnIndex, :) = [0, maximum];
end

fig = figure('Color', 'w', 'Units', 'inches', ...
    'Position', [0.2, 0.2, 16, 14]);
% Match Figure 4E's physical inter-panel spacing.
columnLeft = [0.125, 0.5486];
columnWidth = 0.275;
rowBottom = [0.575, 0.115];
rowHeight = 0.3143;
axesHandles = gobjects(2, 2);

for rowIndex = 1:2
    for columnIndex = 1:2
        ax = axes(fig, 'Position', [columnLeft(columnIndex), ...
            rowBottom(rowIndex), columnWidth, rowHeight]);
        axesHandles(rowIndex, columnIndex) = ax;
        values = loaded.singularMaps.(rowFields{rowIndex}). ...
            (columnFields{columnIndex});
        imagesc(ax, values);
        axis(ax, 'image');
        set(ax, 'YDir', 'normal', 'FontName', 'Arial', 'FontSize', 26, ...
            'LineWidth', 1.4, 'Box', 'on', 'TickDir', 'out', ...
            'XTick', [1 10 20 30 40], 'YTick', [1 10 20 30 40]);
        clim(ax, columnLimits(columnIndex, :));
        colormap(ax, jet(256));
        ax.Toolbar.Visible = 'off';
        hold(ax, 'on');
        for boundary = [10.5 20.5 30.5]
            plot(ax, [boundary boundary], [0.5 40.5], '-', ...
                'Color', [0.45 0.45 0.45], 'LineWidth', 0.7);
            plot(ax, [0.5 40.5], [boundary boundary], '-', ...
                'Color', [0.45 0.45 0.45], 'LineWidth', 0.7);
        end
        hold(ax, 'off');
        if rowIndex == 1
            title(ax, columnTitles{columnIndex}, ...
                'FontName', 'Arial', 'FontSize', 26, 'FontWeight', 'bold');
            ax.XTickLabel = [];
        else
            xlabel(ax, 'Map column (pixel)', 'FontSize', 26);
        end
        if columnIndex == 1
            ylabel(ax, jacobianLabels{rowIndex}, 'FontSize', 30, ...
                'FontWeight', 'bold', 'Interpreter', 'tex');
        else
            ax.YTickLabel = [];
        end
    end
end

for columnIndex = 1:2
    ax = axesHandles(1, columnIndex);
    bar = colorbar(ax);
    ax.Position = [columnLeft(columnIndex), rowBottom(1), ...
        columnWidth, rowHeight];
    bar.Units = 'normalized';
    bar.Position = [columnLeft(columnIndex) + columnWidth + 0.010, ...
        rowBottom(1), 0.016, rowHeight];
    bar.FontName = 'Arial';
    bar.FontSize = 26;
    bar.LineWidth = 1.2;
    bar.Title.String = 'sp/s';
    bar.Title.FontSize = 26;
    bar.Title.FontWeight = 'bold';
end

savefig(fig, outputFig);
exportgraphics(fig, outputPdf, 'ContentType', 'vector', ...
    'BackgroundColor', 'white');
close(fig);
fprintf('Saved %s.\n', outputPdf);
end
