function audit = verify_ji_plus3_atlas_population( ...
        sourceRoot, stagingRoot, populationName)
% Verify one staged J_I atlas while keeping only one FIG open at a time.

stem = sprintf( ...
    '%s_population_J_I_16condition_full_eigenspectrum_maps.fig', ...
    populationName);
sourceFigure = openfig(fullfile(sourceRoot, stem), 'invisible');
[sourceData, sourceLimits] = local_ordered_maps(sourceFigure);
close(sourceFigure);

stagedFigure = openfig(fullfile(stagingRoot, stem), 'invisible');
[stagedData, stagedLimits] = local_ordered_maps(stagedFigure);
close(stagedFigure);

maximumDifferences = zeros(64, 1);
for mapIndex = 1:64
    maximumDifferences(mapIndex) = max(abs( ...
        sourceData{mapIndex}(:) - stagedData{mapIndex}(:)));
end
changedIndices = find(maximumDifferences > 1e-12).';
if ~isequal(changedIndices, [43 59])
    error('Plus3AtlasAudit:ChangedPanels', ...
        '%s changed map indices are [%s], expected [43 59].', ...
        populationName, num2str(changedIndices));
end
maximumLimitDifference = max(abs(sourceLimits(:) - stagedLimits(:)));
if maximumLimitDifference >= 1e-12
    error('Plus3AtlasAudit:ColorLimits', ...
        '%s color limits changed by %.6g.', ...
        populationName, maximumLimitDifference);
end

audit = table(string(populationName), maximumDifferences(43), ...
    maximumDifferences(59), maximumLimitDifference, ...
    'VariableNames', {'population', 'angle15Contrast66Difference', ...
    'angle22p5Contrast66Difference', 'maximumColorLimitDifference'});
fprintf([ ...
    '%s: only indices [43 59] changed, differences %.6g/%.6g; ' ...
    'color limits identical.\n'], populationName, ...
    maximumDifferences(43), maximumDifferences(59));
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
    error('Plus3AtlasAudit:MapCount', ...
        'Expected 64 heatmaps, found %d.', numel(mapAxes));
end

data = cell(64, 1);
limits = zeros(64, 2);
for axisIndex = 1:64
    imageHandle = findobj(mapAxes(axisIndex), 'Type', 'Image');
    data{axisIndex} = double(imageHandle.CData);
    limits(axisIndex, :) = double(mapAxes(axisIndex).CLim);
end
end
