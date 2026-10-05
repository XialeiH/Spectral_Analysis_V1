function result = run_ji_literal_cluster_tolerance_task( ...
        taskIndex, existingDataRoot, outputRoot, tolerance)
% Compute the literal J_I top cluster, including its structural zero modes.

if nargin < 1 || isempty(taskIndex)
    taskIndex = str2double(getenv('JI_CLUSTER_TASK_INDEX'));
    if ~isfinite(taskIndex)
        taskIndex = str2double(getenv('SLURM_ARRAY_TASK_ID'));
    end
end
if nargin < 2 || isempty(existingDataRoot)
    existingDataRoot = getenv('JI_CLUSTER_EXISTING_DATA_ROOT');
end
if nargin < 3 || isempty(outputRoot)
    outputRoot = getenv('JI_CLUSTER_OUTPUT_ROOT');
end
if nargin < 4 || isempty(tolerance)
    tolerance = 5e-2;
end

conditions = [15.0 66; 22.5 66];
if ~isscalar(taskIndex) || ~isfinite(taskIndex) || ...
        taskIndex ~= round(taskIndex) || taskIndex < 1 || ...
        taskIndex > size(conditions, 1)
    error('JICluster:TaskIndex', 'Task index must be 1 or 2.');
end
if ~isfolder(existingDataRoot)
    error('JICluster:DataRoot', ...
        'Missing existing data root: %s', existingDataRoot);
end
if isempty(outputRoot)
    error('JICluster:OutputRoot', 'An output root is required.');
end
if ~exist(outputRoot, 'dir'); mkdir(outputRoot); end

angleDeg = conditions(taskIndex, 1);
contrast = conditions(taskIndex, 2);
angleTag = local_angle_tag(angleDeg);
existingFile = fullfile(existingDataRoot, sprintf( ...
    'full_eigenspectrum_maps_angle%s_contrast%d.mat', ...
    angleTag, contrast));
if ~isfile(existingFile)
    error('JICluster:ExistingFile', ...
        'Missing compact condition file: %s', existingFile);
end

existing = load(existingFile, 'conditionData');
conditionData = existing.conditionData;
sourceFile = conditionData.SourceFiles.DynamicL6;
cacheFile = conditionData.SourceFiles.EigenvalueCache;
source = load(sourceFile, 'Section4');
cache = load(cacheFile, 'eigenvalueData');

jFull = sparse(source.Section4.A);
dimension = size(jFull, 1);
populationSize = dimension / 3;
mapSide = round(sqrt(populationSize));
if ~isequal(size(jFull), [4800 4800]) || mapSide^2 ~= populationSize
    error('JICluster:Dimension', 'Expected a 4800-dimensional Jacobian.');
end
iRows = 2 * populationSize + (1:populationSize);
activeMatrix = full(jFull(iRows, iRows));
cachedEigenvalues = cache.eigenvalueData.J_I(:);

topReal = max(real(cachedEigenvalues));
threshold = topReal - tolerance;
selectionPadding = 100 * eps(max(abs(topReal), 1));
zeroTolerance = 1e-12;
structuralZeroCount = sum(abs(cachedEigenvalues) <= zeroTolerance);
literalSelectedCount = sum(real(cachedEigenvalues) >= ...
    threshold - selectionPadding);

[activeVectors, activeValues] = eig(activeMatrix, 'vector');
[~, order] = sortrows([real(activeValues), imag(activeValues)], [-1 -2]);
activeValues = activeValues(order);
activeVectors = activeVectors(:, order);
selected = find(real(activeValues) >= threshold - selectionPadding);
activeSelectedCount = numel(selected);
if structuralZeroCount + activeSelectedCount ~= literalSelectedCount
    error('JICluster:CountAudit', ...
        ['Structural zeros (%d) plus active modes (%d) do not equal ' ...
         'the cached literal selection count (%d).'], ...
        structuralZeroCount, activeSelectedCount, literalSelectedCount);
end

selectedVectors = activeVectors(:, selected);
selectedVectors = selectedVectors ./ max(vecnorm(selectedVectors), eps);
[activeBasis, triangularFactor] = qr(selectedVectors, 0);
rankTolerance = max(size(triangularFactor)) * ...
    eps(max(norm(triangularFactor, 2), 1));
activeRank = sum(abs(diag(triangularFactor)) > rankTolerance);
activeBasis = activeBasis(:, 1:activeRank);

iEnvelope = reshape(sqrt(sum(abs(activeBasis).^2, 2)), ...
    mapSide, mapSide);
wC = conditionData.CWeight;
wS = 1 - wC;
maps = struct();
maps.S = ones(mapSide);
maps.C = ones(mapSide);
maps.I = iEnvelope;
maps.E = sqrt(wS^2 + wC^2) * ones(mapSide);

selectedValues = activeValues(selected);
residuals = vecnorm(activeMatrix * selectedVectors - ...
    selectedVectors .* reshape(selectedValues, 1, [])) ./ ...
    max((norm(activeMatrix, 'fro') + abs(selectedValues.')) .* ...
    vecnorm(selectedVectors), eps);

result = struct();
result.AngleDeg = angleDeg;
result.Contrast = contrast;
result.Tolerance = tolerance;
result.TopRealEigenvalue = topReal;
result.SelectionThreshold = threshold;
result.StructuralZeroCount = structuralZeroCount;
result.ActiveSelectedCount = activeSelectedCount;
result.LiteralSelectedCount = literalSelectedCount;
result.ActiveRank = activeRank;
result.LiteralClusterRank = structuralZeroCount + activeRank;
result.CurrentTolerance = conditionData.ClusterTolerance;
result.CurrentClusterCount = ...
    conditionData.Jacobians.J_I.TopClusterCount;
result.SelectedActiveEigenvalues = selectedValues;
result.EnvelopeMaps = maps;
result.MaximumActiveResidual = max(residuals);
result.SourceConditionFile = existingFile;
result.SourceGeometryFile = sourceFile;
result.SourceEigenvalueCacheFile = cacheFile;

outputFile = fullfile(outputRoot, sprintf( ...
    'J_I_literal_cluster_tol5e-2_angle%s_contrast%d.mat', ...
    angleTag, contrast));
save(outputFile, 'result', '-v7');
fprintf(['Saved literal J_I cluster: angle %.2f, contrast %d, ' ...
    'top %.8g, threshold %.8g, modes %d (%d zeros + %d active), ' ...
    'rank %d -> %s\n'], angleDeg, contrast, topReal, threshold, ...
    literalSelectedCount, structuralZeroCount, activeSelectedCount, ...
    result.LiteralClusterRank, outputFile);
end

function tag = local_angle_tag(angleDeg)
tag = strrep(sprintf('%.2f', angleDeg), '.', 'p');
end
