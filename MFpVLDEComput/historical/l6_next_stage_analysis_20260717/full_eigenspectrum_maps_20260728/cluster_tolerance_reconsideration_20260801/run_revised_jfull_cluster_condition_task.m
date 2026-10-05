function conditionData = run_revised_jfull_cluster_condition_task( ...
        taskIndex, existingDataRoot, outputRoot)
% Replace only the J_full top-cluster fields using the settled condition rules.

if nargin < 1 || isempty(taskIndex)
    taskIndex = str2double(getenv('REVISED_CLUSTER_TASK_INDEX'));
    if ~isfinite(taskIndex)
        taskIndex = str2double(getenv('SLURM_ARRAY_TASK_ID'));
    end
end
if nargin < 2 || isempty(existingDataRoot)
    existingDataRoot = getenv('REVISED_CLUSTER_EXISTING_DATA_ROOT');
end
if nargin < 3 || isempty(outputRoot)
    outputRoot = getenv('REVISED_CLUSTER_OUTPUT_ROOT');
end

angles = [0 7.5 15 22.5];
contrasts = [19 42 66 100];
conditionCount = numel(angles) * numel(contrasts);
if ~isscalar(taskIndex) || ~isfinite(taskIndex) || taskIndex ~= round(taskIndex) || ...
        taskIndex < 1 || taskIndex > conditionCount
    error('RevisedCluster:TaskIndex', 'Task index must be an integer from 1 to 16.');
end
if ~isfolder(existingDataRoot)
    error('RevisedCluster:DataRoot', 'Missing existing data root: %s', existingDataRoot);
end
if isempty(outputRoot)
    error('RevisedCluster:OutputRoot', 'An output root is required.');
end
if ~exist(outputRoot, 'dir'); mkdir(outputRoot); end

[contrastIndex, angleIndex] = ind2sub( ...
    [numel(contrasts), numel(angles)], taskIndex);
angleDeg = angles(angleIndex);
contrast = contrasts(contrastIndex);
angleTag = local_angle_tag(angleDeg);
fileName = sprintf('full_eigenspectrum_maps_angle%s_contrast%d.mat', ...
    angleTag, contrast);
inputFile = fullfile(existingDataRoot, fileName);
outputFile = fullfile(outputRoot, fileName);
if ~isfile(inputFile)
    error('RevisedCluster:InputFile', 'Missing condition file: %s', inputFile);
end

rule = local_selection_rule(angleDeg, contrast);
if strcmp(rule.Type, 'original')
    copyfile(inputFile, outputFile, 'f');
    loaded = load(inputFile, 'conditionData');
    conditionData = loaded.conditionData;
    fprintf('Copied unchanged original cluster: angle %.2f, contrast %d.\n', ...
        angleDeg, contrast);
    return
end

loaded = load(inputFile, 'conditionData');
conditionData = loaded.conditionData;
sourceFile = conditionData.SourceFiles.DynamicL6;
cacheFile = conditionData.SourceFiles.EigenvalueCache;
if ~isfile(sourceFile) || ~isfile(cacheFile)
    error('RevisedCluster:SourceFiles', 'A required source file is missing.');
end
source = load(sourceFile, 'Section4');
cache = load(cacheFile, 'eigenvalueData');
jFull = sparse(source.Section4.A);
cachedEigenvalues = cache.eigenvalueData.J_baseline(:);
dimension = size(jFull, 1);
if ~isequal(size(jFull), [4800 4800]) || numel(cachedEigenvalues) ~= dimension
    error('RevisedCluster:Dimension', 'Expected a 4800-dimensional J_full spectrum.');
end

topReal = max(real(cachedEigenvalues));
selectionPadding = 100 * eps(max(abs(topReal), 1));
if strcmp(rule.Type, 'tolerance')
    selectedCount = sum(real(cachedEigenvalues) >= ...
        topReal - (rule.Value + selectionPadding));
else
    selectedCount = rule.Value;
end
eigsCount = min(dimension - 2, max(16, selectedCount + 8));
options = struct('tol', 1e-10, 'maxit', 5000, ...
    'p', min(dimension, max(80, 2*eigsCount + 16)), ...
    'disp', 0, 'isreal', true);
[vectors, diagonal, flag] = eigs(jFull, eigsCount, 'largestreal', options);
if flag ~= 0
    error('RevisedCluster:EigsConvergence', 'eigs returned flag %d.', flag);
end
values = diag(diagonal);
[~, order] = sortrows([real(values), imag(values)], [-1 -2]);
values = values(order);
vectors = vectors(:, order);
if abs(real(values(1)) - topReal) > 1e-8
    error('RevisedCluster:TopEigenvalueAudit', ...
        'Recovered top real eigenvalue differs from the cached spectrum.');
end
if strcmp(rule.Type, 'tolerance')
    selected = find(real(values) >= ...
        topReal - (rule.Value + selectionPadding), selectedCount, 'first');
else
    selected = 1:selectedCount;
end
if numel(selected) ~= selectedCount
    error('RevisedCluster:SelectionIncomplete', ...
        'Recovered %d of %d requested modes.', numel(selected), selectedCount);
end

selectedValues = values(selected);
selectedVectors = vectors(:, selected);
selectedVectors = selectedVectors ./ max(vecnorm(selectedVectors), eps);
[basis, triangularFactor] = qr(selectedVectors, 0);
rankTolerance = max(size(triangularFactor)) * ...
    eps(max(norm(triangularFactor, 2), 1));
clusterRank = sum(abs(diag(triangularFactor)) > rankTolerance);
basis = basis(:, 1:clusterRank);

populationSize = dimension / 3;
mapSide = round(sqrt(populationSize));
wC = conditionData.CWeight;
clusterMaps = local_cluster_envelope_maps( ...
    basis, populationSize, mapSide, wC);
residuals = vecnorm(jFull * selectedVectors - ...
    selectedVectors .* reshape(selectedValues, 1, [])) ./ ...
    max((norm(jFull, 'fro') + abs(selectedValues.')) .* ...
    vecnorm(selectedVectors), eps);
if max(residuals) > 1e-10
    error('RevisedCluster:Residual', ...
        'Selected eigenvector residual %.3g exceeds tolerance.', max(residuals));
end

analysis = conditionData.Jacobians.J_full;
analysis.TopClusterEnvelopeMaps = clusterMaps;
analysis.TopClusterEigenvalues = selectedValues;
analysis.TopClusterCount = selectedCount;
analysis.TopClusterRank = clusterRank;
conditionData.Jacobians.J_full = analysis;
save(outputFile, 'conditionData', '-v7');
fprintf(['Revised only J_full cluster: angle %.2f, contrast %d, %s %.4g, ' ...
    'count %d, rank %d -> %s\n'], angleDeg, contrast, rule.Type, ...
    rule.Value, selectedCount, clusterRank, outputFile);
end

function rule = local_selection_rule(angleDeg, contrast)
if angleDeg == 7.5 && ismember(contrast, [42 66 100])
    rule = struct('Type', 'tolerance', 'Value', 5e-3);
elseif angleDeg == 15
    if contrast == 66
        modeCount = 10;
    else
        modeCount = 6;
    end
    rule = struct('Type', 'fixed_count', 'Value', modeCount);
elseif angleDeg == 22.5
    if ismember(contrast, [19 100])
        modeCount = 6;
    else
        modeCount = 5;
    end
    rule = struct('Type', 'fixed_count', 'Value', modeCount);
else
    rule = struct('Type', 'original', 'Value', 1e-3);
end
end

function maps = local_cluster_envelope_maps(basis, populationSize, mapSide, wC)
wS = 1 - wC;
s = basis(1:populationSize, :);
c = basis(populationSize + (1:populationSize), :);
i = basis(2*populationSize + (1:populationSize), :);
e = wS*s + wC*c;
maps = struct('S', reshape(sqrt(sum(abs(s).^2, 2)), mapSide, mapSide), ...
    'C', reshape(sqrt(sum(abs(c).^2, 2)), mapSide, mapSide), ...
    'I', reshape(sqrt(sum(abs(i).^2, 2)), mapSide, mapSide), ...
    'E', reshape(sqrt(sum(abs(e).^2, 2)), mapSide, mapSide));
end

function tag = local_angle_tag(angleDeg)
tag = strrep(sprintf('%.2f', angleDeg), '.', 'p');
end
