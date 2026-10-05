function manifest = replot_report_full_eigenspectrum_atlases( ...
        sourceRoot, outputRoot, requestedPopulation, requestedJacobian, ...
        jiClusterOverrideRoot, jiClusterOverrideVariant, ...
        preserveStoredColorLimits)
% Replot the finalized 12 atlases with coordinated report formatting.

if nargin < 1 || isempty(sourceRoot)
    error('ReportAtlas:SourceRoot', 'A source FIG directory is required.');
end
if nargin < 2 || isempty(outputRoot)
    error('ReportAtlas:OutputRoot', 'An output directory is required.');
end
if nargin < 3
    requestedPopulation = '';
end
if nargin < 4
    requestedJacobian = '';
end
if nargin < 5
    jiClusterOverrideRoot = '';
end
if nargin < 6 || isempty(jiClusterOverrideVariant)
    jiClusterOverrideVariant = 'literal5e2';
end
if nargin < 7 || isempty(preserveStoredColorLimits)
    preserveStoredColorLimits = false;
end
if ~isfolder(sourceRoot)
    error('ReportAtlas:MissingSourceRoot', 'Missing source root: %s', sourceRoot);
end
if ~exist(outputRoot, 'dir'); mkdir(outputRoot); end

jacobianNames = {'J_full', 'J_6', 'J_I'};
populationNames = {'E', 'S', 'C', 'I'};
if ~isempty(requestedPopulation) && ...
        ~ismember(requestedPopulation, populationNames)
    error('ReportAtlas:Population', ...
        'Unsupported population: %s', requestedPopulation);
end
if ~isempty(requestedJacobian) && ...
        ~ismember(requestedJacobian, jacobianNames)
    error('ReportAtlas:Jacobian', ...
        'Unsupported Jacobian: %s', requestedJacobian);
end
loadedPopulationNames = populationNames;
if ~isempty(requestedPopulation)
    loadedPopulationNames = {requestedPopulation};
end
angles = repelem([0 7.5 15 22.5], 4);
contrasts = repmat([19 42 66 100], 1, 4);

atlases = struct();
for populationIndex = 1:numel(loadedPopulationNames)
    populationName = loadedPopulationNames{populationIndex};
    for jacobianIndex = 1:numel(jacobianNames)
        jacobianName = jacobianNames{jacobianIndex};
        sourceFile = fullfile(sourceRoot, sprintf( ...
            '%s_population_%s_16condition_full_eigenspectrum_maps.fig', ...
            populationName, jacobianName));
        atlases.(populationName).(jacobianName) = ...
            local_extract_atlas(sourceFile);
    end
end

sharedLimits = struct();
modeColorPercentile = 99;
for populationIndex = 1:numel(loadedPopulationNames)
    populationName = loadedPopulationNames{populationIndex};
    limits = zeros(4, 2);
    for columnIndex = 1:4
        values = [];
        for jacobianIndex = 1:numel(jacobianNames)
            atlas = atlases.(populationName).(jacobianNames{jacobianIndex});
            for rowIndex = 1:numel(atlas.Rows)
                values = [values; atlas.Rows{rowIndex}.MapData{columnIndex}(:)]; %#ok<AGROW>
            end
        end
        values = values(isfinite(values));
        if columnIndex == 1
            limits(columnIndex, :) = [min(values) max(values)];
        elseif ismember(columnIndex, [2 4])
            absoluteLimit = prctile(abs(values), modeColorPercentile);
            limits(columnIndex, :) = [-absoluteLimit absoluteLimit];
        else
            limits(columnIndex, :) = ...
                [0 prctile(values, modeColorPercentile)];
        end
        limits(columnIndex, :) = local_expand_equal_limits( ...
            limits(columnIndex, :));
    end
    sharedLimits.(populationName) = limits;
end

if preserveStoredColorLimits || ~isempty(jiClusterOverrideRoot)
    for populationIndex = 1:numel(loadedPopulationNames)
        populationName = loadedPopulationNames{populationIndex};
        sourceAtlas = atlases.(populationName).J_full;
        storedLimits = sourceAtlas.Rows{1}.MapLimits;
        sharedLimits.(populationName) = storedLimits;
        for jacobianIndex = 1:numel(jacobianNames)
            sourceAtlas = atlases.(populationName).( ...
                jacobianNames{jacobianIndex});
            for rowIndex = 1:numel(sourceAtlas.Rows)
                if max(abs(sourceAtlas.Rows{rowIndex}.MapLimits(:) - ...
                        storedLimits(:))) >= 1e-12
                    error('ReportAtlas:StoredColorLimits', ...
                        ['%s source atlases do not share identical map ' ...
                         'limits across Jacobians and rows.'], ...
                        populationName);
                end
            end
        end
    end
end
if ~isempty(jiClusterOverrideRoot)
    atlases = local_apply_ji_cluster_overrides(atlases, ...
        loadedPopulationNames, angles, contrasts, jiClusterOverrideRoot, ...
        jiClusterOverrideVariant);
end

plottedJacobianNames = jacobianNames;
if ~isempty(requestedJacobian)
    plottedJacobianNames = {requestedJacobian};
end
manifestRows = cell(numel(plottedJacobianNames) * ...
    numel(loadedPopulationNames), 11);
manifestIndex = 0;
for jacobianIndex = 1:numel(plottedJacobianNames)
    jacobianName = plottedJacobianNames{jacobianIndex};
    for populationIndex = 1:numel(loadedPopulationNames)
        populationName = loadedPopulationNames{populationIndex};
        atlas = atlases.(populationName).(jacobianName);
        limits = sharedLimits.(populationName);
        outputs = local_plot_atlas(atlas, populationName, jacobianName, ...
            limits, angles, contrasts, outputRoot);
        manifestIndex = manifestIndex + 1;
        manifestRows(manifestIndex, :) = {string(jacobianName), ...
            string(populationName), string(outputs.FigureFile), ...
            string(outputs.PdfFile), limits(1, 1), limits(1, 2), ...
            max(abs(limits(2, :))), limits(3, 2), ...
            max(abs(limits(4, :))), modeColorPercentile, 600};
    end
end

manifest = cell2table(manifestRows, 'VariableNames', { ...
    'jacobian', 'population', 'figureFile', 'pdfFile', ...
    'sharedFixedPointMinimum', 'sharedFixedPointMaximum', ...
    'sharedTopModeAbsoluteLimit', 'sharedClusterEnvelopeMaximum', ...
    'sharedSingularOutputAbsoluteLimit', 'modeColorPercentile', ...
    'exportResolutionDpi'});
manifestName = 'report_atlas_manifest.tsv';
if ~isempty(requestedPopulation) && ~isempty(requestedJacobian)
    manifestName = sprintf('report_atlas_manifest_%s_%s.tsv', ...
        requestedPopulation, requestedJacobian);
end
writetable(manifest, fullfile(outputRoot, manifestName), ...
    'FileType', 'text', 'Delimiter', '\t');
fprintf('Saved %d coordinated report atlases to %s.\n', ...
    height(manifest), outputRoot);
end

function atlases = local_apply_ji_cluster_overrides(atlases, ...
        populationNames, angles, contrasts, overrideRoot, overrideVariant)
if ~isfolder(overrideRoot)
    error('ReportAtlas:OverrideRoot', ...
        'Missing J_I cluster override root: %s', overrideRoot);
end

overrideConditions = [15 66; 22.5 66];
for conditionIndex = 1:size(overrideConditions, 1)
    angleDeg = overrideConditions(conditionIndex, 1);
    contrast = overrideConditions(conditionIndex, 2);
    rowIndex = find(angles == angleDeg & contrasts == contrast);
    if numel(rowIndex) ~= 1
        error('ReportAtlas:OverrideRow', ...
            'Could not identify one atlas row for angle %.2f, contrast %d.', ...
            angleDeg, contrast);
    end

    angleTag = strrep(sprintf('%.2f', angleDeg), '.', 'p');
    switch overrideVariant
        case 'literal5e2'
            overrideFile = fullfile(overrideRoot, sprintf( ...
                'J_I_literal_cluster_tol5e-2_angle%s_contrast%d.mat', ...
                angleTag, contrast));
            mapField = 'EnvelopeMaps';
        case 'plus3'
            overrideFile = fullfile(overrideRoot, sprintf( ...
                'J_I_expanded_cluster_plus3_angle%s_contrast%d.mat', ...
                angleTag, contrast));
            mapField = 'ExpandedEnvelopeMaps';
        otherwise
            error('ReportAtlas:OverrideVariant', ...
                'Unsupported J_I cluster override variant: %s', ...
                overrideVariant);
    end
    if ~isfile(overrideFile)
        error('ReportAtlas:OverrideFile', ...
            'Missing J_I cluster override: %s', overrideFile);
    end
    loaded = load(overrideFile, 'result');
    if strcmp(overrideVariant, 'plus3')
        expectedCount = 17 - (angleDeg == 22.5);
        if loaded.result.ExpandedCount ~= expectedCount
            error('ReportAtlas:OverrideModeCount', ...
                'Expected %d expanded modes in %s, found %d.', ...
                expectedCount, overrideFile, loaded.result.ExpandedCount);
        end
    end
    for populationIndex = 1:numel(populationNames)
        populationName = populationNames{populationIndex};
        mapContainer = loaded.result.(mapField);
        map = double(mapContainer.(populationName));
        if ~isequal(size(map), [40 40])
            error('ReportAtlas:OverrideMap', ...
                'Expected a 40-by-40 %s override map in %s.', ...
                populationName, overrideFile);
        end
        atlases.(populationName).J_I.Rows{rowIndex}.MapData{3} = map;
    end
end
end

function atlas = local_extract_atlas(sourceFile)
if ~isfile(sourceFile)
    error('ReportAtlas:MissingSource', 'Missing source FIG: %s', sourceFile);
end

sourceFigure = openfig(sourceFile, 'invisible');
cleanup = onCleanup(@() close(sourceFigure));
allAxes = findall(sourceFigure, 'Type', 'axes');
isPanel = false(size(allAxes));
for axisIndex = 1:numel(allAxes)
    isPanel(axisIndex) = ...
        ~isempty(findobj(allAxes(axisIndex), 'Type', 'Image')) || ...
        ~isempty(findobj(allAxes(axisIndex), 'Type', 'Scatter')) || ...
        ~isempty(findobj(allAxes(axisIndex), 'Type', 'Bar'));
end
panelAxes = allAxes(isPanel);
if numel(panelAxes) ~= 16 * 6
    error('ReportAtlas:PanelCount', ...
        'Expected 96 panel axes in %s, found %d.', sourceFile, numel(panelAxes));
end

pixelPositions = zeros(numel(panelAxes), 4);
for axisIndex = 1:numel(panelAxes)
    pixelPositions(axisIndex, :) = getpixelposition(panelAxes(axisIndex), true);
end
centerY = pixelPositions(:, 2) + pixelPositions(:, 4) / 2;
[~, rowOrder] = sort(centerY, 'descend');
panelAxes = panelAxes(rowOrder);

rows = cell(16, 1);
for rowIndex = 1:16
    rowAxes = panelAxes((rowIndex - 1) * 6 + (1:6));
    xPositions = arrayfun(@(axisHandle) ...
        getpixelposition(axisHandle, true), rowAxes, 'UniformOutput', false);
    xPositions = cellfun(@(position) position(1), xPositions);
    [~, columnOrder] = sort(xPositions);
    rowAxes = rowAxes(columnOrder);
    rows{rowIndex} = local_extract_row(rowAxes, sourceFile, rowIndex);
end
atlas = struct('Rows', {rows}, 'SourceFile', sourceFile);
end

function row = local_extract_row(rowAxes, sourceFile, rowIndex)
row = struct();
row.MapData = cell(4, 1);
row.MapLimits = zeros(4, 2);
for columnIndex = 1:4
    imageHandle = findobj(rowAxes(columnIndex), 'Type', 'Image');
    if numel(imageHandle) ~= 1
        error('ReportAtlas:Image', ...
            'Expected one heatmap in row %d, column %d of %s.', ...
            rowIndex, columnIndex, sourceFile);
    end
    row.MapData{columnIndex} = double(imageHandle.CData);
    row.MapLimits(columnIndex, :) = double(rowAxes(columnIndex).CLim);
end

scatterHandle = findobj(rowAxes(5), 'Type', 'Scatter');
barHandle = findobj(rowAxes(6), 'Type', 'Bar');
if numel(scatterHandle) ~= 1 || numel(barHandle) ~= 1
    error('ReportAtlas:Spectrum', ...
        'Could not identify spectrum and histogram in row %d of %s.', ...
        rowIndex, sourceFile);
end
row.SpectrumX = double(scatterHandle.XData(:));
row.SpectrumY = double(scatterHandle.YData(:));
row.SpectrumXLimits = double(rowAxes(5).XLim);
row.SpectrumYLimits = double(rowAxes(5).YLim);
row.HistogramX = double(barHandle.XData(:));
row.HistogramY = double(barHandle.YData(:));
row.HistogramXLimits = double(rowAxes(6).XLim);
row.HistogramYLimits = double(rowAxes(6).YLim);
end

function outputs = local_plot_atlas(atlas, populationName, jacobianName, ...
        sharedMapLimits, angles, contrasts, outputRoot)
rowCount = numel(atlas.Rows);
fig = figure('Visible', 'off', 'Color', 'w', 'Units', 'pixels', ...
    'Position', [100 100 1800 2600]);
colormap(fig, jet(256));

mapStarts = [0.065 0.190 0.315 0.440];
mapWidth = 0.078;
colorbarGap = 0.003;
colorbarWidth = 0.006;
spectrumStart = 0.615;
histogramStart = 0.755;
spectralWidth = 0.105;
topAxes = 0.930;
bottomAxes = 0.022;
rowStep = (topAxes - bottomAxes) / rowCount;
panelHeight = 0.0448;

columnTitles = {'Fixed Point', 'Top Eigenvector', 'Top Eigencluster', ...
    'Top Singular Vector', 'Eigenspectrum', 'Eigenvalue Histogram'};
columnCenters = [mapStarts + mapWidth/2, ...
    spectrumStart + spectralWidth/2, histogramStart + spectralWidth/2];
for columnIndex = 1:numel(columnTitles)
    annotation(fig, 'textbox', ...
        [columnCenters(columnIndex) - 0.07 0.956 0.14 0.018], ...
        'String', columnTitles{columnIndex}, 'EdgeColor', 'none', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
        'FontSize', 12, 'FontWeight', 'bold', 'Interpreter', 'none', ...
        'FitBoxToText', 'off', 'Margin', 0);
end

for rowIndex = 1:rowCount
    row = atlas.Rows{rowIndex};
    rowY = topAxes - rowIndex * rowStep + ...
        (rowStep - panelHeight) / 2;
    coordinates = local_sampling_coordinates(angles(rowIndex));
    for columnIndex = 1:4
        axisHandle = axes(fig, 'Position', ...
            [mapStarts(columnIndex) rowY mapWidth panelHeight]);
        map = row.MapData{columnIndex};
        imagesc(axisHandle, map);
        axis(axisHandle, 'image');
        axis(axisHandle, 'off');
        clim(axisHandle, sharedMapLimits(columnIndex, :));
        local_add_hypercolumn_grid(axisHandle, size(map));

        sampleValues = [map(coordinates(1, 1), coordinates(1, 2)), ...
            local_five_pixel_surround_mean(map, coordinates(2, :)), ...
            map(coordinates(3, 1), coordinates(3, 2))];
        titleHandle = title(axisHandle, local_value_triplet(sampleValues), ...
            'Interpreter', 'none', 'FontSize', 6.8, ...
            'FontWeight', 'normal');
        titleHandle.Position(2) = 1.015;

        colorbarHandle = colorbar(axisHandle, 'Location', 'eastoutside');
        axisHandle.Position = ...
            [mapStarts(columnIndex) rowY mapWidth panelHeight];
        colorbarHandle.Position = [mapStarts(columnIndex) + mapWidth + ...
            colorbarGap rowY colorbarWidth panelHeight];
        local_format_colorbar(colorbarHandle, ...
            sharedMapLimits(columnIndex, :));

        if columnIndex == 1
            text(axisHandle, -0.27, 0.5, sprintf( ...
                'angle %.1f deg\ncontrast %d', ...
                angles(rowIndex), contrasts(rowIndex)), ...
                'Units', 'normalized', 'HorizontalAlignment', 'right', ...
                'VerticalAlignment', 'middle', 'FontSize', 8.2, ...
                'FontWeight', 'bold', 'Interpreter', 'none', ...
                'Clipping', 'off');
        end
    end

    spectrumAxis = axes(fig, 'Position', ...
        [spectrumStart rowY spectralWidth panelHeight]);
    local_spectrum_panel(spectrumAxis, row);

    histogramAxis = axes(fig, 'Position', ...
        [histogramStart rowY spectralWidth panelHeight]);
    local_histogram_panel(histogramAxis, row);
end

annotation(fig, 'textbox', [0.05 0.972 0.86 0.018], ...
    'String', sprintf([ ...
        '%s population; %s; ' ...
        '(x, y, z) = (Pref, Surround 5-pixel mean, Ortho)'], ...
        populationName, jacobianName), ...
    'Interpreter', 'none', 'FontSize', 16, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
    'EdgeColor', 'none', 'Margin', 0);

stem = sprintf('%s_population_%s_16condition_full_eigenspectrum_maps', ...
    populationName, jacobianName);
figureFile = fullfile(outputRoot, [stem '.fig']);
pdfFile = fullfile(outputRoot, [stem '.pdf']);
savefig(fig, figureFile);
exportgraphics(fig, pdfFile, 'ContentType', 'vector', 'Resolution', 600);
close(fig);
pdfInfo = dir(pdfFile);
if isempty(pdfInfo) || pdfInfo.bytes < 1e5
    error('ReportAtlas:BlankPdf', ...
        ['Vector export produced an empty or truncated PDF for %s. ' ...
         'Re-export its saved FIG in a fresh MATLAB process.'], stem);
end

outputs = struct('FigureFile', figureFile, 'PdfFile', pdfFile);
fprintf('Saved coordinated %s population; %s atlas.\n', ...
    populationName, jacobianName);
end

function value = local_five_pixel_surround_mean(map, center)
row = center(1);
column = center(2);
neighbors = [row column; row - 1 column; row + 1 column; ...
    row column - 1; row column + 1];
if any(neighbors(:) < 1) || any(neighbors(:, 1) > size(map, 1)) || ...
        any(neighbors(:, 2) > size(map, 2))
    error('ReportAtlas:SurroundNeighborhood', ...
        'The five-pixel surround neighborhood exceeds the map boundary.');
end
indices = sub2ind(size(map), neighbors(:, 1), neighbors(:, 2));
value = mean(map(indices));
end

function coordinates = local_sampling_coordinates(angleDeg)
if angleDeg == 0
    coordinates = [5 10; 5 7; 5 1];
elseif angleDeg == 7.5
    coordinates = [4 10; 5 7; 6 1];
elseif angleDeg == 15
    coordinates = [3 10; 5 5; 10 1];
elseif angleDeg == 22.5
    coordinates = [1 10; 5 5; 10 1];
else
    error('ReportAtlas:Angle', 'Unsupported angle %.6g.', angleDeg);
end
end

function label = local_value_triplet(values)
parts = arrayfun(@local_format_value, values, 'UniformOutput', false);
label = sprintf('(%s, %s, %s)', parts{:});
end

function label = local_format_value(value)
if value == 0
    label = '0';
elseif abs(value) >= 1e-3 && abs(value) < 1e3
    label = sprintf('%.3g', value);
else
    label = sprintf('%.2e', value);
end
end

function local_add_hypercolumn_grid(axisHandle, mapSize)
hold(axisHandle, 'on');
for boundary = 10.5:10:(mapSize(2) - 0.5)
    line(axisHandle, [boundary boundary], [0.5 mapSize(1) + 0.5], ...
        'Color', [0.45 0.45 0.45], 'LineWidth', 0.35, ...
        'HitTest', 'off', 'HandleVisibility', 'off');
end
for boundary = 10.5:10:(mapSize(1) - 0.5)
    line(axisHandle, [0.5 mapSize(2) + 0.5], [boundary boundary], ...
        'Color', [0.45 0.45 0.45], 'LineWidth', 0.35, ...
        'HitTest', 'off', 'HandleVisibility', 'off');
end
hold(axisHandle, 'off');
end

function local_format_colorbar(colorbarHandle, limits)
colorbarHandle.FontSize = 5.8;
colorbarHandle.LineWidth = 0.45;
colorbarHandle.TickDirection = 'out';
colorbarHandle.Ticks = linspace(limits(1), limits(2), 3);
colorbarHandle.TickLabels = arrayfun(@local_format_value, ...
    colorbarHandle.Ticks, 'UniformOutput', false);
end

function local_spectrum_panel(axisHandle, row)
scatter(axisHandle, row.SpectrumX, row.SpectrumY, 4, ...
    [0 0.4470 0.7410], 'filled', 'MarkerFaceAlpha', 0.58, ...
    'MarkerEdgeAlpha', 0.20);
hold(axisHandle, 'on');
xline(axisHandle, 0, ':', 'Color', [0.3 0.3 0.3], 'LineWidth', 0.5);
if row.SpectrumXLimits(1) <= 1 && row.SpectrumXLimits(2) >= 1
    xline(axisHandle, 1, '--', 'Color', [0.75 0.15 0.15], ...
        'LineWidth', 0.5);
end
hold(axisHandle, 'off');
xlim(axisHandle, row.SpectrumXLimits);
ylim(axisHandle, row.SpectrumYLimits);
pbaspect(axisHandle, [1 1 1]);
local_format_spectral_axis(axisHandle, ...
    'Re(\lambda)', 'Im(\lambda)', row.SpectrumXLimits, ...
    row.SpectrumYLimits);
end

function local_histogram_panel(axisHandle, row)
bar(axisHandle, row.HistogramX, row.HistogramY, 1, ...
    'FaceColor', [0 0.4470 0.7410], 'EdgeColor', 'none');
hold(axisHandle, 'on');
xline(axisHandle, 0, ':', 'Color', [0.3 0.3 0.3], 'LineWidth', 0.5);
if row.HistogramXLimits(1) <= 1 && row.HistogramXLimits(2) >= 1
    xline(axisHandle, 1, '--', 'Color', [0.75 0.15 0.15], ...
        'LineWidth', 0.5);
end
hold(axisHandle, 'off');
xlim(axisHandle, row.HistogramXLimits);
ylim(axisHandle, row.HistogramYLimits);
pbaspect(axisHandle, [1 1 1]);
local_format_spectral_axis(axisHandle, ...
    'Re(\lambda)', 'Count', row.HistogramXLimits, ...
    row.HistogramYLimits);
end

function local_format_spectral_axis(axisHandle, xLabel, yLabel, ...
        xLimits, yLimits)
set(axisHandle, 'Box', 'on', 'LineWidth', 0.55, ...
    'FontSize', 5.8, 'Color', 'w', 'TickDir', 'out');
axisHandle.XTick = linspace(xLimits(1), xLimits(2), 3);
axisHandle.YTick = linspace(yLimits(1), yLimits(2), 3);
axisHandle.XTickLabel = arrayfun(@local_format_value, ...
    axisHandle.XTick, 'UniformOutput', false);
axisHandle.YTickLabel = arrayfun(@local_format_value, ...
    axisHandle.YTick, 'UniformOutput', false);
xlabel(axisHandle, xLabel, 'Interpreter', 'tex', 'FontSize', 6.2);
ylabel(axisHandle, yLabel, 'Interpreter', 'tex', 'FontSize', 6.2);
end

function limits = local_expand_equal_limits(limits)
if limits(2) <= limits(1)
    scale = max(abs(limits(1)), 1);
    limits = limits(1) + [-1 1] * scale * eps;
end
end
