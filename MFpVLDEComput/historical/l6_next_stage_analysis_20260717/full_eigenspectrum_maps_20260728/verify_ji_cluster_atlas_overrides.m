function audit = verify_ji_cluster_atlas_overrides( ...
        sourceRoot, stagingRoot, requestedPopulation)
% Verify that only the two requested J_I cluster panels changed.

populationNames = {'E', 'S', 'C', 'I'};
if nargin >= 3 && ~isempty(requestedPopulation)
    if ~ismember(requestedPopulation, populationNames)
        error('ReportAtlasAudit:Population', ...
            'Unsupported population: %s', requestedPopulation);
    end
    populationNames = {requestedPopulation};
end
expectedChangedIndices = [43 59];
auditRows = cell(numel(populationNames), 5);
for populationIndex = 1:numel(populationNames)
    populationName = populationNames{populationIndex};
    stem = sprintf( ...
        '%s_population_J_I_16condition_full_eigenspectrum_maps.fig', ...
        populationName);
    sourceFigure = openfig(fullfile(sourceRoot, stem), 'invisible');
    sourceCleanup = onCleanup(@() close(sourceFigure));
    stagedFigure = openfig(fullfile(stagingRoot, stem), 'invisible');
    stagedCleanup = onCleanup(@() close(stagedFigure));

    [sourceData, sourceLimits] = local_ordered_maps(sourceFigure);
    [stagedData, stagedLimits] = local_ordered_maps(stagedFigure);
    maximumDifferences = zeros(64, 1);
    for mapIndex = 1:64
        maximumDifferences(mapIndex) = max(abs( ...
            sourceData{mapIndex}(:) - stagedData{mapIndex}(:)));
    end
    changedIndices = find(maximumDifferences > 1e-12).';
    if ~isequal(changedIndices, expectedChangedIndices)
        error('ReportAtlasAudit:ChangedPanels', ...
            '%s changed map indices are [%s], expected [43 59].', ...
            populationName, num2str(changedIndices));
    end
    maximumLimitDifference = max(abs( ...
        sourceLimits(:) - stagedLimits(:)));
    if maximumLimitDifference >= 1e-12
        error('ReportAtlasAudit:ColorLimits', ...
            '%s color limits changed by %.6g.', ...
            populationName, maximumLimitDifference);
    end
    auditRows(populationIndex, :) = {populationName, ...
        changedIndices, maximumDifferences(43), ...
        maximumDifferences(59), maximumLimitDifference};
    fprintf([ ...
        '%s: changed indices [43 59], differences %.6g/%.6g; ' ...
        'color limits identical.\n'], populationName, ...
        maximumDifferences(43), maximumDifferences(59));
    clear stagedCleanup sourceCleanup
end

audit = cell2table(auditRows, 'VariableNames', { ...
    'population', 'changedMapIndices', 'angle15Contrast66Difference', ...
    'angle22p5Contrast66Difference', 'maximumColorLimitDifference'});
end

function [data, limits] = local_ordered_maps(fig)
allAxes = findall(fig, 'Type', 'axes');
isMap = false(size(allAxes));
for axisIndex = 1:numel(allAxes)
    isMap(axisIndex) = ...
        isscalar(findobj(allAxes(axisIndex), 'Type', 'Image'));
end
mapAxes = allAxes(isMap);
positions = zeros(numel(mapAxes), 2);
for axisIndex = 1:numel(mapAxes)
    position = getpixelposition(mapAxes(axisIndex), true);
    positions(axisIndex, :) = ...
        [-(position(2) + position(4) / 2), position(1)];
end
[~, order] = sortrows(positions, [1 2]);
mapAxes = mapAxes(order);
if numel(mapAxes) ~= 64
    error('ReportAtlasAudit:MapCount', ...
        'Expected 64 heatmaps, found %d.', numel(mapAxes));
end

data = cell(numel(mapAxes), 1);
limits = zeros(numel(mapAxes), 2);
for axisIndex = 1:numel(mapAxes)
    imageHandle = findobj(mapAxes(axisIndex), 'Type', 'Image');
    data{axisIndex} = double(imageHandle.CData);
    limits(axisIndex, :) = double(mapAxes(axisIndex).CLim);
end
end
