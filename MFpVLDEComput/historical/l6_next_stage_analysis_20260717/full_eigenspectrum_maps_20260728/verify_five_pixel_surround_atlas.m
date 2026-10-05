function audit = verify_five_pixel_surround_atlas( ...
        sourceRoot, stagingRoot, populationName, jacobianName)
% Verify that replotting changed labels but not numerical panel data.

stem = sprintf('%s_population_%s_16condition_full_eigenspectrum_maps.fig', ...
    populationName, jacobianName);
sourceFigure = openfig(fullfile(sourceRoot, stem), 'invisible');
source = local_extract(sourceFigure);
close(sourceFigure);
stagedFigure = openfig(fullfile(stagingRoot, stem), 'invisible');
staged = local_extract(stagedFigure);
close(stagedFigure);

maximumMapDifference = 0;
maximumColorLimitDifference = 0;
maximumSpectrumDifference = 0;
maximumHistogramDifference = 0;
angles = repelem([0 7.5 15 22.5], 4);
for rowIndex = 1:16
    coordinates = local_sampling_coordinates(angles(rowIndex));
    for columnIndex = 1:4
        maximumMapDifference = max(maximumMapDifference, max(abs( ...
            source.Rows{rowIndex}.MapData{columnIndex}(:) - ...
            staged.Rows{rowIndex}.MapData{columnIndex}(:))));
        maximumColorLimitDifference = max(maximumColorLimitDifference, ...
            max(abs(source.Rows{rowIndex}.MapLimits(columnIndex, :) - ...
            staged.Rows{rowIndex}.MapLimits(columnIndex, :))));

        map = staged.Rows{rowIndex}.MapData{columnIndex};
        values = [map(coordinates(1, 1), coordinates(1, 2)), ...
            local_five_pixel_surround_mean(map, coordinates(2, :)), ...
            map(coordinates(3, 1), coordinates(3, 2))];
        expectedTitle = local_value_triplet(values);
        if ~strcmp(staged.Rows{rowIndex}.MapTitles{columnIndex}, expectedTitle)
            error('FivePixelAudit:MapTitle', ...
                ['%s %s row %d column %d title is "%s"; ' ...
                 'expected "%s".'], populationName, jacobianName, ...
                rowIndex, columnIndex, ...
                staged.Rows{rowIndex}.MapTitles{columnIndex}, expectedTitle);
        end
    end

    maximumSpectrumDifference = max(maximumSpectrumDifference, max([ ...
        max(abs(source.Rows{rowIndex}.SpectrumX - ...
            staged.Rows{rowIndex}.SpectrumX)), ...
        max(abs(source.Rows{rowIndex}.SpectrumY - ...
            staged.Rows{rowIndex}.SpectrumY)), ...
        max(abs(source.Rows{rowIndex}.SpectrumXLimits - ...
            staged.Rows{rowIndex}.SpectrumXLimits)), ...
        max(abs(source.Rows{rowIndex}.SpectrumYLimits - ...
            staged.Rows{rowIndex}.SpectrumYLimits))]));
    maximumHistogramDifference = max(maximumHistogramDifference, max([ ...
        max(abs(source.Rows{rowIndex}.HistogramX - ...
            staged.Rows{rowIndex}.HistogramX)), ...
        max(abs(source.Rows{rowIndex}.HistogramY - ...
            staged.Rows{rowIndex}.HistogramY)), ...
        max(abs(source.Rows{rowIndex}.HistogramXLimits - ...
            staged.Rows{rowIndex}.HistogramXLimits)), ...
        max(abs(source.Rows{rowIndex}.HistogramYLimits - ...
            staged.Rows{rowIndex}.HistogramYLimits))]));
end

tolerance = 1e-12;
if any([maximumMapDifference, maximumColorLimitDifference, ...
        maximumSpectrumDifference, maximumHistogramDifference] >= tolerance)
    error('FivePixelAudit:NumericalDifference', ...
        ['%s %s changed numerical panel data: map %.3g, CLim %.3g, ' ...
         'spectrum %.3g, histogram %.3g.'], ...
        populationName, jacobianName, maximumMapDifference, ...
        maximumColorLimitDifference, maximumSpectrumDifference, ...
        maximumHistogramDifference);
end

audit = table(string(populationName), string(jacobianName), ...
    maximumMapDifference, maximumColorLimitDifference, ...
    maximumSpectrumDifference, maximumHistogramDifference, ...
    'VariableNames', {'population', 'jacobian', 'maximumMapDifference', ...
    'maximumColorLimitDifference', 'maximumSpectrumDifference', ...
    'maximumHistogramDifference'});
fprintf('%s %s: numerical panels identical; all 64 titles verified.\n', ...
    populationName, jacobianName);
end

function atlas = local_extract(fig)
allAxes = findall(fig, 'Type', 'axes');
isPanel = false(size(allAxes));
for axisIndex = 1:numel(allAxes)
    isPanel(axisIndex) = ...
        ~isempty(findobj(allAxes(axisIndex), 'Type', 'Image')) || ...
        ~isempty(findobj(allAxes(axisIndex), 'Type', 'Scatter')) || ...
        ~isempty(findobj(allAxes(axisIndex), 'Type', 'Bar'));
end
panelAxes = allAxes(isPanel);
positions = zeros(numel(panelAxes), 4);
for axisIndex = 1:numel(panelAxes)
    positions(axisIndex, :) = ...
        getpixelposition(panelAxes(axisIndex), true);
end
centerY = positions(:, 2) + positions(:, 4) / 2;
[~, rowOrder] = sort(centerY, 'descend');
panelAxes = panelAxes(rowOrder);
if numel(panelAxes) ~= 96
    error('FivePixelAudit:PanelCount', ...
        'Expected 96 panel axes, found %d.', numel(panelAxes));
end

rows = cell(16, 1);
for rowIndex = 1:16
    rowAxes = panelAxes((rowIndex - 1) * 6 + (1:6));
    xPositions = zeros(6, 1);
    for columnIndex = 1:6
        position = getpixelposition(rowAxes(columnIndex), true);
        xPositions(columnIndex) = position(1);
    end
    [~, columnOrder] = sort(xPositions);
    rowAxes = rowAxes(columnOrder);
    row = struct('MapData', {cell(4, 1)}, ...
        'MapTitles', {cell(4, 1)}, 'MapLimits', zeros(4, 2));
    for columnIndex = 1:4
        imageHandle = findobj(rowAxes(columnIndex), 'Type', 'Image');
        row.MapData{columnIndex} = double(imageHandle.CData);
        row.MapTitles{columnIndex} = char(rowAxes(columnIndex).Title.String);
        row.MapLimits(columnIndex, :) = double(rowAxes(columnIndex).CLim);
    end
    scatterHandle = findobj(rowAxes(5), 'Type', 'Scatter');
    row.SpectrumX = double(scatterHandle.XData(:));
    row.SpectrumY = double(scatterHandle.YData(:));
    row.SpectrumXLimits = double(rowAxes(5).XLim);
    row.SpectrumYLimits = double(rowAxes(5).YLim);
    barHandle = findobj(rowAxes(6), 'Type', 'Bar');
    row.HistogramX = double(barHandle.XData(:));
    row.HistogramY = double(barHandle.YData(:));
    row.HistogramXLimits = double(rowAxes(6).XLim);
    row.HistogramYLimits = double(rowAxes(6).YLim);
    rows{rowIndex} = row;
end
atlas = struct('Rows', {rows});
end

function value = local_five_pixel_surround_mean(map, center)
row = center(1);
column = center(2);
indices = sub2ind(size(map), ...
    [row; row - 1; row + 1; row; row], ...
    [column; column; column; column - 1; column + 1]);
value = mean(map(indices));
end

function coordinates = local_sampling_coordinates(angleDeg)
if angleDeg == 0
    coordinates = [5 10; 5 7; 5 1];
elseif angleDeg == 7.5
    coordinates = [4 10; 5 7; 6 1];
elseif angleDeg == 15
    coordinates = [3 10; 5 5; 10 1];
else
    coordinates = [1 10; 5 5; 10 1];
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
