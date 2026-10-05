function result = run_jfull_c_cluster_tolerance_task(taskIndex, existingDataRoot, outputRoot)
% Recompute selected J_full C-population eigencluster envelopes at two tolerances.

if nargin < 1 || isempty(taskIndex)
    taskIndex = str2double(getenv('CLUSTER_TASK_INDEX'));
    if ~isfinite(taskIndex)
        taskIndex = str2double(getenv('SLURM_ARRAY_TASK_ID'));
    end
end
if nargin < 2 || isempty(existingDataRoot)
    existingDataRoot = getenv('CLUSTER_EXISTING_DATA_ROOT');
end
if nargin < 3 || isempty(outputRoot)
    outputRoot = getenv('CLUSTER_OUTPUT_ROOT');
end

profile = getenv('CLUSTER_PROFILE');
if isempty(profile); profile = 'initial'; end
[conditions, tolerances] = local_profile(profile);

if ~isscalar(taskIndex) || ~isfinite(taskIndex) || taskIndex ~= round(taskIndex) || ...
        taskIndex < 1 || taskIndex > size(conditions, 1)
    error('ClusterTolerance:TaskIndex', ...
        'Task index must be an integer from 1 to %d.', size(conditions, 1));
end
if ~isfolder(existingDataRoot)
    error('ClusterTolerance:DataRoot', 'Missing existing data root: %s', existingDataRoot);
end
if isempty(outputRoot)
    error('ClusterTolerance:OutputRoot', 'An output root is required.');
end
if ~exist(outputRoot, 'dir'); mkdir(outputRoot); end

angleDeg = conditions(taskIndex, 1);
contrast = conditions(taskIndex, 2);
angleTag = local_angle_tag(angleDeg);
existingFile = fullfile(existingDataRoot, sprintf( ...
    'full_eigenspectrum_maps_angle%s_contrast%d.mat', angleTag, contrast));
if ~isfile(existingFile)
    error('ClusterTolerance:ExistingFile', 'Missing compact condition file: %s', existingFile);
end

existing = load(existingFile, 'conditionData');
sourceFile = existing.conditionData.SourceFiles.DynamicL6;
cacheFile = existing.conditionData.SourceFiles.EigenvalueCache;
if ~isfile(sourceFile)
    error('ClusterTolerance:SourceFile', 'Missing dynamic-L6 geometry file: %s', sourceFile);
end
if ~isfile(cacheFile)
    error('ClusterTolerance:CacheFile', 'Missing eigenvalue cache: %s', cacheFile);
end

source = load(sourceFile, 'Section4');
cache = load(cacheFile, 'eigenvalueData');
jFull = sparse(source.Section4.A);
eigenvalues = cache.eigenvalueData.J_baseline(:);
dimension = size(jFull, 1);
if ~isequal(size(jFull), [4800 4800]) || numel(eigenvalues) ~= dimension
    error('ClusterTolerance:Dimension', 'Expected a 4800-dimensional J_full spectrum.');
end

topReal = max(real(eigenvalues));
selectionPadding = 100 * eps(max(abs(topReal), 1));
counts = arrayfun(@(tol) sum(real(eigenvalues) >= ...
    topReal - (tol + selectionPadding)), tolerances);
maximumExpectedCount = max(counts);
requestedCount = min(dimension - 2, max(16, maximumExpectedCount + 8));
maximumCount = min(dimension - 2, max(192, maximumExpectedCount + 32));

while true
    options = struct('tol', 1e-10, 'maxit', 5000, ...
        'p', min(dimension, max(80, 2*requestedCount + 16)), ...
        'disp', 0, 'isreal', true);
    [vectors, diagonal, flag] = eigs(jFull, requestedCount, 'largestreal', options);
    if flag ~= 0
        error('ClusterTolerance:EigsConvergence', 'eigs returned flag %d.', flag);
    end
    values = diag(diagonal);
    [~, order] = sortrows([real(values), imag(values)], [-1 -2]);
    values = values(order);
    vectors = vectors(:, order);
    recovered = sum(real(values) >= ...
        topReal - (max(tolerances) + selectionPadding));
    if recovered >= maximumExpectedCount || requestedCount >= maximumCount
        break
    end
    requestedCount = min(maximumCount, ...
        max(2*requestedCount, maximumExpectedCount + 8));
end
if recovered < maximumExpectedCount
    error('ClusterTolerance:ClusterIncomplete', ...
        'Recovered %d of %d modes at tolerance %.4g.', ...
        recovered, maximumExpectedCount, max(tolerances));
end

populationSize = dimension / 3;
mapSide = round(sqrt(populationSize));
result = struct();
result.AngleDeg = angleDeg;
result.Contrast = contrast;
result.Tolerances = tolerances;
result.Profile = profile;
result.TopRealEigenvalue = topReal;
result.SourceFile = sourceFile;
result.EigenvalueCacheFile = cacheFile;
result.Cluster = repmat(struct(), 1, numel(tolerances));

for toleranceIndex = 1:numel(tolerances)
    tolerance = tolerances(toleranceIndex);
    expectedCount = counts(toleranceIndex);
    selected = find(real(values) >= ...
        topReal - (tolerance + selectionPadding), expectedCount, 'first');
    if numel(selected) ~= expectedCount
        error('ClusterTolerance:SelectionIncomplete', ...
            'Recovered %d of %d modes at tolerance %.4g.', ...
            numel(selected), expectedCount, tolerance);
    end

    selectedVectors = vectors(:, selected);
    selectedVectors = selectedVectors ./ max(vecnorm(selectedVectors), eps);
    [basis, triangularFactor] = qr(selectedVectors, 0);
    rankTolerance = max(size(triangularFactor)) * ...
        eps(max(norm(triangularFactor, 2), 1));
    clusterRank = sum(abs(diag(triangularFactor)) > rankTolerance);
    basis = basis(:, 1:clusterRank);
    cBasis = basis(populationSize + (1:populationSize), :);
    cMap = reshape(sqrt(sum(abs(cBasis).^2, 2)), mapSide, mapSide);

    selectedValues = values(selected);
    residuals = vecnorm(jFull * selectedVectors - ...
        selectedVectors .* reshape(selectedValues, 1, [])) ./ ...
        max((norm(jFull, 'fro') + abs(selectedValues.')) .* ...
        vecnorm(selectedVectors), eps);

    result.Cluster(toleranceIndex).Tolerance = tolerance;
    result.Cluster(toleranceIndex).Count = expectedCount;
    result.Cluster(toleranceIndex).Rank = clusterRank;
    result.Cluster(toleranceIndex).Eigenvalues = selectedValues;
    result.Cluster(toleranceIndex).CEnvelopeMap = cMap;
    result.Cluster(toleranceIndex).MaximumResidual = max(residuals);
end

outputFile = fullfile(outputRoot, sprintf( ...
    'Jfull_C_cluster_tolerance_angle%s_contrast%d.mat', angleTag, contrast));
save(outputFile, 'result', '-v7');
fprintf(['Saved J_full C-cluster tolerance comparison: angle %.2f, contrast %d, ' ...
    'counts %d/%d -> %s\n'], angleDeg, contrast, counts(1), counts(2), outputFile);
end

function tag = local_angle_tag(angleDeg)
tag = strrep(sprintf('%.2f', angleDeg), '.', 'p');
end

function [conditions, tolerances] = local_profile(profile)
switch profile
    case 'initial'
        conditions = [ ...
             7.5,  42;  7.5,  66;  7.5, 100; ...
            15.0,  19; 15.0,  42; 15.0,  66; 15.0, 100; ...
            22.5,  19; 22.5,  42; 22.5,  66; 22.5, 100];
        tolerances = [5e-4, 5e-3];
    case 'remaining'
        conditions = [ ...
            15.0,  19; 15.0,  42; 15.0,  66; 15.0, 100; ...
            22.5,  19; 22.5,  42; 22.5,  66; 22.5, 100];
        tolerances = [1e-4, 1e-2];
    otherwise
        error('ClusterTolerance:Profile', 'Unsupported profile: %s', profile);
end
end
