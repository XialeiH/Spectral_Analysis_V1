function [summary, conditions] = ...
        compute_fixedpoint_mode_similarity_from_atlases( ...
        atlasRoot, outputRoot)
% Compute fixed-point alignment with cluster envelopes and singular modes.

if ~isfolder(atlasRoot)
    error('FixedModeSimilarity:AtlasRoot', ...
        'Missing atlas root: %s', atlasRoot);
end
if ~exist(outputRoot, 'dir'); mkdir(outputRoot); end

jacobianNames = {'J_full', 'J_6', 'J_I'};
populationNames = {'S', 'C', 'I'};
angles = repelem([0 7.5 15 22.5], 4).';
contrasts = repmat([19 42 66 100], 1, 4).';

atlas = struct();
for populationIndex = 1:numel(populationNames)
    populationName = populationNames{populationIndex};
    for jacobianIndex = 1:numel(jacobianNames)
        jacobianName = jacobianNames{jacobianIndex};
        sourceFile = fullfile(atlasRoot, sprintf( ...
            '%s_population_%s_16condition_full_eigenspectrum_maps.fig', ...
            populationName, jacobianName));
        atlas.(populationName).(jacobianName) = ...
            local_extract_mode_maps(sourceFile);
    end
end

clusterSimilarity = zeros(16, 3);
singularSimilarity = zeros(16, 3);
maximumFixedPointMismatch = 0;
for conditionIndex = 1:16
    fixedPoint = local_full_state_vector( ...
        atlas, populationNames, 'J_full', 'FixedPoint', conditionIndex);
    for jacobianIndex = 1:3
        jacobianName = jacobianNames{jacobianIndex};
        comparisonFixedPoint = local_full_state_vector( ...
            atlas, populationNames, jacobianName, ...
            'FixedPoint', conditionIndex);
        maximumFixedPointMismatch = max(maximumFixedPointMismatch, ...
            max(abs(fixedPoint - comparisonFixedPoint)));

        clusterEnvelope = local_full_state_vector( ...
            atlas, populationNames, jacobianName, ...
            'ClusterEnvelope', conditionIndex);
        singularMode = local_full_state_vector( ...
            atlas, populationNames, jacobianName, ...
            'SingularMode', conditionIndex);
        clusterSimilarity(conditionIndex, jacobianIndex) = ...
            local_cosine(fixedPoint, clusterEnvelope, false);
        singularSimilarity(conditionIndex, jacobianIndex) = ...
            local_cosine(fixedPoint, singularMode, true);
    end
end
if maximumFixedPointMismatch >= 1e-12
    error('FixedModeSimilarity:FixedPointMismatch', ...
        'Fixed-point maps differ across Jacobian atlases by %.6g.', ...
        maximumFixedPointMismatch);
end

conditions = table(angles, contrasts, ...
    clusterSimilarity(:, 1), clusterSimilarity(:, 2), ...
    clusterSimilarity(:, 3), singularSimilarity(:, 1), ...
    singularSimilarity(:, 2), singularSimilarity(:, 3), ...
    'VariableNames', {'angleDeg', 'contrast', ...
    'fixedClusterJFull', 'fixedClusterJ6', 'fixedClusterJI', ...
    'fixedSingularJFull', 'fixedSingularJ6', 'fixedSingularJI'});

summaryRows = cell(6, 6);
summaryIndex = 0;
for metricIndex = 1:2
    if metricIndex == 1
        values = clusterSimilarity;
        quantity = 'Fixed point vs top eigencluster';
    else
        values = singularSimilarity;
        quantity = 'Fixed point vs top singular output';
    end
    for jacobianIndex = 1:3
        summaryIndex = summaryIndex + 1;
        column = values(:, jacobianIndex);
        summaryRows(summaryIndex, :) = {quantity, ...
            jacobianNames{jacobianIndex}, mean(column), median(column), ...
            min(column), max(column)};
    end
end
summary = cell2table(summaryRows, 'VariableNames', { ...
    'quantity', 'jacobian', 'mean', 'median', 'minimum', 'maximum'});

writetable(summary, fullfile(outputRoot, ...
    'fixedpoint_mode_similarity_summary.tsv'), ...
    'FileType', 'text', 'Delimiter', '\t');
writetable(conditions, fullfile(outputRoot, ...
    'fixedpoint_mode_similarity_conditions.tsv'), ...
    'FileType', 'text', 'Delimiter', '\t');
save(fullfile(outputRoot, 'fixedpoint_mode_similarity_results.mat'), ...
    'summary', 'conditions', 'maximumFixedPointMismatch', '-v7');
disp(summary);
end

function maps = local_extract_mode_maps(sourceFile)
if ~isfile(sourceFile)
    error('FixedModeSimilarity:SourceFile', ...
        'Missing atlas FIG: %s', sourceFile);
end
fig = openfig(sourceFile, 'invisible');
cleanup = onCleanup(@() close(fig));
allAxes = findall(fig, 'Type', 'axes');
isPanel = false(size(allAxes));
for axisIndex = 1:numel(allAxes)
    isPanel(axisIndex) = ...
        ~isempty(findobj(allAxes(axisIndex), 'Type', 'Image')) || ...
        ~isempty(findobj(allAxes(axisIndex), 'Type', 'Scatter')) || ...
        ~isempty(findobj(allAxes(axisIndex), 'Type', 'Bar'));
end
panelAxes = allAxes(isPanel);
if numel(panelAxes) ~= 96
    error('FixedModeSimilarity:PanelCount', ...
        'Expected 96 panels in %s, found %d.', ...
        sourceFile, numel(panelAxes));
end

pixelPositions = zeros(numel(panelAxes), 4);
for axisIndex = 1:numel(panelAxes)
    pixelPositions(axisIndex, :) = ...
        getpixelposition(panelAxes(axisIndex), true);
end
centerY = pixelPositions(:, 2) + pixelPositions(:, 4) / 2;
[~, rowOrder] = sort(centerY, 'descend');
panelAxes = panelAxes(rowOrder);

maps = struct('FixedPoint', {cell(16, 1)}, ...
    'ClusterEnvelope', {cell(16, 1)}, ...
    'SingularMode', {cell(16, 1)});
for rowIndex = 1:16
    rowAxes = panelAxes((rowIndex - 1) * 6 + (1:6));
    xPositions = zeros(6, 1);
    for axisIndex = 1:6
        position = getpixelposition(rowAxes(axisIndex), true);
        xPositions(axisIndex) = position(1);
    end
    [~, columnOrder] = sort(xPositions);
    rowAxes = rowAxes(columnOrder);
    maps.FixedPoint{rowIndex} = ...
        local_image_data(rowAxes(1), sourceFile, rowIndex, 1);
    maps.ClusterEnvelope{rowIndex} = ...
        local_image_data(rowAxes(3), sourceFile, rowIndex, 3);
    maps.SingularMode{rowIndex} = ...
        local_image_data(rowAxes(4), sourceFile, rowIndex, 4);
end
end

function data = local_image_data(axisHandle, sourceFile, rowIndex, columnIndex)
imageHandle = findobj(axisHandle, 'Type', 'Image');
if ~isscalar(imageHandle)
    error('FixedModeSimilarity:Image', ...
        'Expected one image in row %d, column %d of %s.', ...
        rowIndex, columnIndex, sourceFile);
end
data = double(imageHandle.CData);
end

function vector = local_full_state_vector( ...
        atlas, populationNames, jacobianName, mapName, conditionIndex)
parts = cell(numel(populationNames), 1);
for populationIndex = 1:numel(populationNames)
    populationName = populationNames{populationIndex};
    mapSet = atlas.(populationName).(jacobianName).(mapName);
    parts{populationIndex} = mapSet{conditionIndex}(:);
end
vector = vertcat(parts{:});
end

function value = local_cosine(x, y, phaseInvariant)
denominator = norm(x) * norm(y);
if denominator == 0
    error('FixedModeSimilarity:ZeroVector', ...
        'Cannot compute cosine similarity with a zero vector.');
end
innerProduct = x' * y;
if phaseInvariant
    innerProduct = abs(innerProduct);
end
value = real(innerProduct / denominator);
end
