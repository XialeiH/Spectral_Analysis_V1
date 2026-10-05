function outputs = make_2_22_1_shared_log_colorbar(sourceFigure, outputRoot)
% Preserve Figure 2.22 and change only to shared column-wise log colorbars.

if nargin < 1 || ~isfile(sourceFigure)
    error('SharedLogColorbar:Source', 'A valid source FIG is required.');
end
if nargin < 2 || isempty(outputRoot)
    error('SharedLogColorbar:OutputRoot', 'An output directory is required.');
end
if ~exist(outputRoot, 'dir'); mkdir(outputRoot); end

fig = openfig(sourceFigure, 'invisible');
dataAxes = findall(fig, 'Type', 'axes');
dataAxes = dataAxes(arrayfun(@(axisHandle) ...
    ~isempty(findall(axisHandle, 'Type', 'image')), dataAxes));
if numel(dataAxes) ~= 12
    close(fig);
    error('SharedLogColorbar:Axes', 'Expected 12 heatmap axes.');
end

positions = vertcat(dataAxes.Position);
[~, columnOrder] = sort(positions(:, 1));
columnIndex = zeros(12, 1);
for rankIndex = 1:12
    columnIndex(columnOrder(rankIndex)) = ceil(rankIndex / 3);
end

columnMaximum = zeros(1, 4);
for axisIndex = 1:12
    imageHandle = findall(dataAxes(axisIndex), 'Type', 'image');
    data = imageHandle(1).CData;
    column = columnIndex(axisIndex);
    columnMaximum(column) = max(columnMaximum(column), max(abs(data(:))));
end
columnMaximum = max(columnMaximum, eps);
linearThreshold = columnMaximum * 1e-3;

colorbars = findall(fig, 'Type', 'ColorBar');
if numel(colorbars) ~= 12
    close(fig);
    error('SharedLogColorbar:Colorbars', 'Expected 12 panel colorbars.');
end

for axisIndex = 1:12
    axisHandle = dataAxes(axisIndex);
    imageHandle = findall(axisHandle, 'Type', 'image');
    data = imageHandle(1).CData;
    column = columnIndex(axisIndex);
    if column <= 2
        displayData = log10(1 + max(data, 0) ./ ...
            linearThreshold(column));
        axisHandle.CLim = [0 log10(1 + columnMaximum(column) ./ ...
            linearThreshold(column))];
    else
        displayData = sign(data) .* log10(1 + abs(data) ./ ...
            linearThreshold(column));
        displayLimit = log10(1 + columnMaximum(column) ./ ...
            linearThreshold(column));
        axisHandle.CLim = [-displayLimit displayLimit];
    end
    imageHandle(1).CData = displayData;

end

colorbarPositions = vertcat(colorbars.Position);
[~, colorbarOrder] = sort(colorbarPositions(:, 1));
colorbarColumn = zeros(12, 1);
for rankIndex = 1:12
    colorbarColumn(colorbarOrder(rankIndex)) = ceil(rankIndex / 3);
end
for colorbarIndex = 1:12
    column = colorbarColumn(colorbarIndex);
    if column <= 2
        originalTicks = columnMaximum(column) * [0 0.001 0.01 0.1 1];
        displayTicks = log10(1 + originalTicks ./ ...
            linearThreshold(column));
    else
        originalTicks = columnMaximum(column) * ...
            [-1 -0.1 -0.01 0 0.01 0.1 1];
        displayTicks = sign(originalTicks) .* log10( ...
            1 + abs(originalTicks) ./ linearThreshold(column));
    end
    colorbars(colorbarIndex).Ticks = displayTicks;
    colorbars(colorbarIndex).TickLabels = compose('%.2g', originalTicks);
end

layout = findall(fig, 'Type', 'tiledlayout');
if isscalar(layout)
    oldTitle = string(layout.Title.String);
    oldTitle(1) = regexprep(oldTitle(1), '^2\.19:', '2.22.1:');
    oldTitle(2) = ['Shared logarithmic color scale within each column; ' ...
        'all maps and mode definitions unchanged'];
    title(layout, cellstr(oldTitle), 'FontSize', 15, ...
        'FontWeight', 'bold', 'Interpreter', 'none');
end

stem = ['2.22.1_Normalized_Left_Right_Cluster_and_Singular_Maps_' ...
    'SharedLogColor_Angle0_Contrast100'];
outputs = struct();
outputs.FigureFile = fullfile(outputRoot, [stem '.fig']);
outputs.PdfFile = fullfile(outputRoot, [stem '.pdf']);
savefig(fig, outputs.FigureFile);
exportgraphics(fig, outputs.PdfFile, ...
    'ContentType', 'vector', 'Resolution', 600);
close(fig);
fprintf('Saved shared-log colorbar figure to %s.\n', outputs.PdfFile);
end
