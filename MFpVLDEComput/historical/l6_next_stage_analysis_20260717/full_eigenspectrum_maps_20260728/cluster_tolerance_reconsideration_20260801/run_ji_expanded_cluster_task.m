function result = run_ji_expanded_cluster_task( ...
        taskIndex, existingDataRoot, outputRoot, extraModeCount)
% Expand the original J_I cluster by an exact number of leading modes.

if nargin < 1 || isempty(taskIndex)
    taskIndex = str2double(getenv('SLURM_ARRAY_TASK_ID'));
end
if nargin < 2 || isempty(existingDataRoot)
    existingDataRoot = getenv('JI_CLUSTER_EXISTING_DATA_ROOT');
end
if nargin < 3 || isempty(outputRoot)
    outputRoot = getenv('JI_CLUSTER_OUTPUT_ROOT');
end
if nargin < 4 || isempty(extraModeCount)
    extraModeCount = 5;
end

conditions = [15 66; 22.5 66];
if ~isscalar(taskIndex) || taskIndex ~= round(taskIndex) || ...
        taskIndex < 1 || taskIndex > size(conditions, 1)
    error('JIExpandedCluster:TaskIndex', 'Task index must be 1 or 2.');
end
if ~isfolder(existingDataRoot)
    error('JIExpandedCluster:DataRoot', ...
        'Missing existing data root: %s', existingDataRoot);
end
if ~exist(outputRoot, 'dir'); mkdir(outputRoot); end

angleDeg = conditions(taskIndex, 1);
contrast = conditions(taskIndex, 2);
angleTag = strrep(sprintf('%.2f', angleDeg), '.', 'p');
existingFile = fullfile(existingDataRoot, sprintf( ...
    'full_eigenspectrum_maps_angle%s_contrast%d.mat', ...
    angleTag, contrast));
loaded = load(existingFile, 'conditionData');
conditionData = loaded.conditionData;
source = load(conditionData.SourceFiles.DynamicL6, 'Section4');

jFull = sparse(source.Section4.A);
dimension = size(jFull, 1);
n = dimension / 3;
side = round(sqrt(n));
eRows = 1:(2*n);
iRows = 2*n + (1:n);
jI = sparse(dimension, dimension);
jI(:, iRows) = jFull(:, iRows);
activeMatrix = full(jI(iRows, iRows));

originalCount = conditionData.Jacobians.J_I.TopClusterCount;
expandedCount = originalCount + extraModeCount;
[activeVectors, activeValues] = eig(activeMatrix, 'vector');
[~, order] = sortrows([real(activeValues), imag(activeValues)], [-1 -2]);
activeValues = activeValues(order);
activeVectors = activeVectors(:, order);
selectedValues = activeValues(1:expandedCount);
if any(abs(selectedValues) <= 1e-12)
    error('JIExpandedCluster:ZeroMode', ...
        'The expanded active cluster unexpectedly reached a zero mode.');
end

lifted = complex(zeros(dimension, expandedCount));
for modeIndex = 1:expandedCount
    vector = activeVectors(:, modeIndex);
    lambda = selectedValues(modeIndex);
    fullVector = [jI(eRows, iRows) * vector / lambda; vector];
    lifted(:, modeIndex) = fullVector / max(norm(fullVector), eps);
end
[basis, triangularFactor] = qr(lifted, 0);
rankTolerance = max(size(triangularFactor)) * ...
    eps(max(norm(triangularFactor, 2), 1));
expandedRank = sum(abs(diag(triangularFactor)) > rankTolerance);
basis = basis(:, 1:expandedRank);

residuals = vecnorm(jI * lifted - ...
    lifted .* reshape(selectedValues, 1, [])) ./ ...
    max((norm(jI, 'fro') + abs(selectedValues.')) .* ...
    vecnorm(lifted), eps);

result = struct();
result.AngleDeg = angleDeg;
result.Contrast = contrast;
result.OriginalCount = originalCount;
result.ExpandedCount = expandedCount;
result.ExtraModeCount = extraModeCount;
result.ExpandedRank = expandedRank;
result.TopRealEigenvalue = real(selectedValues(1));
result.LastIncludedRealEigenvalue = real(selectedValues(end));
result.ImpliedTolerance = result.TopRealEigenvalue - ...
    result.LastIncludedRealEigenvalue;
result.SelectedEigenvalues = selectedValues;
result.OriginalEnvelopeMaps = ...
    conditionData.Jacobians.J_I.TopClusterEnvelopeMaps;
result.ExpandedEnvelopeMaps = local_cluster_envelope_maps( ...
    basis, n, side, conditionData.CWeight);
result.MaximumResidual = max(residuals);
result.SourceConditionFile = existingFile;

outputFile = fullfile(outputRoot, sprintf( ...
    'J_I_expanded_cluster_%s_angle%s_contrast%d.mat', ...
    local_extra_tag(extraModeCount), angleTag, contrast));
save(outputFile, 'result', '-v7');
fprintf([ ...
    'Saved J_I expanded cluster: angle %.1f, C%d, %d -> %d modes, ' ...
    'implied tolerance %.8g, max residual %.3g.\n'], ...
    angleDeg, contrast, originalCount, expandedCount, ...
    result.ImpliedTolerance, result.MaximumResidual);
end

function tag = local_extra_tag(extraModeCount)
if extraModeCount >= 0
    tag = sprintf('plus%d', extraModeCount);
else
    tag = sprintf('minus%d', abs(extraModeCount));
end
end

function maps = local_cluster_envelope_maps(basis, n, side, wC)
wS = 1 - wC;
s = basis(1:n, :);
c = basis(n + (1:n), :);
i = basis(2*n + (1:n), :);
e = wS*s + wC*c;
maps = struct( ...
    'S', reshape(sqrt(sum(abs(s).^2, 2)), side, side), ...
    'C', reshape(sqrt(sum(abs(c).^2, 2)), side, side), ...
    'I', reshape(sqrt(sum(abs(i).^2, 2)), side, side), ...
    'E', reshape(sqrt(sum(abs(e).^2, 2)), side, side));
end
