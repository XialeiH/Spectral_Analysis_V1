function result = compute_normalized_e_heatmap_left_right_task( ...
        taskIndex, sourceRoot, modeSummaryFile, outputRoot)
% Compare unit-L2 E heatmaps for left/right clusters and singular vectors.

if nargin < 1 || isempty(taskIndex)
    taskIndex = str2double(getenv('SLURM_ARRAY_TASK_ID'));
end
if nargin < 2 || isempty(sourceRoot)
    sourceRoot = getenv('NORMALIZED_MAP_SOURCE_ROOT');
end
if nargin < 3 || isempty(modeSummaryFile)
    modeSummaryFile = getenv('NORMALIZED_MAP_MODE_SUMMARY');
end
if nargin < 4 || isempty(outputRoot)
    outputRoot = getenv('NORMALIZED_MAP_OUTPUT_ROOT');
end
if ~isscalar(taskIndex) || taskIndex ~= round(taskIndex) || ...
        taskIndex < 1 || taskIndex > 16
    error('NormalizedMap:TaskIndex', ...
        'Task index must be an integer from 1 to 16.');
end
if ~isfolder(sourceRoot) || ~isfile(modeSummaryFile)
    error('NormalizedMap:Inputs', 'A source root and mode summary are required.');
end
if ~exist(outputRoot, 'dir'); mkdir(outputRoot); end

angles = [0 7.5 15 22.5];
contrasts = [19 42 66 100];
[contrastIndex, angleIndex] = ind2sub([4 4], taskIndex);
angleDeg = angles(angleIndex);
contrast = contrasts(contrastIndex);

dynamicFile = fullfile(sourceRoot, 'geometry_mat', sprintf( ...
    'geometry_sections4_5_h96baseline_L6eqW0p00_contr%d_angle_%.2f.mat', ...
    contrast, angleDeg));
frozenFile = fullfile(sourceRoot, 'geometry_mat', sprintf( ...
    'geometry_sections4_5_h96baseline_L6eqW1p00_contr%d_angle_%.2f.mat', ...
    contrast, angleDeg));
dynamicData = load(dynamicFile, 'Section4');
frozenData = load(frozenFile, 'Section4');
jFull = sparse(dynamicData.Section4.A);
jFrozen = sparse(frozenData.Section4.A);
clear dynamicData frozenData

dimension = size(jFull, 1);
n = dimension / 3;
if dimension ~= 4800 || n ~= round(n)
    error('NormalizedMap:Dimension', 'Expected a 4800-dimensional Jacobian.');
end
eRows = 1:2*n;
iRows = 2*n + (1:n);
j6 = sparse(jFull - jFrozen);
jI = sparse(dimension, dimension);
jI(:, iRows) = jFull(:, iRows);
clear jFrozen

modeSummary = readtable(modeSummaryFile, ...
    'FileType', 'text', 'Delimiter', '\t', 'TextType', 'string');
jacobianNames = {'J_full', 'J_6', 'J_I'};
matrices = {jFull, j6, jI};
rows = cell(3, 12);
conditionMaps = struct();
for jacobianIndex = 1:3
    jacobianName = jacobianNames{jacobianIndex};
    selectedSummary = modeSummary.jacobian == string(jacobianName) & ...
        abs(modeSummary.angleDeg - angleDeg) < 1e-9 & ...
        modeSummary.contrast == contrast;
    if sum(selectedSummary) ~= 1
        error('NormalizedMap:ModeCount', ...
            'Missing unique mode count for %s, angle %.2f, contrast %d.', ...
            jacobianName, angleDeg, contrast);
    end
    modeCount = modeSummary.clusterModeCount(selectedSummary);
    matrix = matrices{jacobianIndex};
    [rightEnvelope, leftEnvelope] = local_cluster_e_maps( ...
        matrix, jacobianName, modeCount, n, eRows, iRows);
    [rightSingular, leftSingular] = local_singular_e_maps(matrix, n);

    [clusterHC, clusterCosine, rightClusterMap, leftClusterMap] = ...
        local_normalized_map_metrics(rightEnvelope, leftEnvelope, false);
    [singularHC, singularCosine, rightSingularMap, leftSingularMap] = ...
        local_normalized_map_metrics(rightSingular, leftSingular, true);
    conditionMaps.(jacobianName) = struct( ...
        'RightCluster', reshape(rightClusterMap, 40, 40), ...
        'LeftCluster', reshape(leftClusterMap, 40, 40), ...
        'RightSingular', reshape(rightSingularMap, 40, 40), ...
        'LeftSingular', reshape(leftSingularMap, 40, 40), ...
        'ClusterModeCount', modeCount);
    rows(jacobianIndex, :) = {string(jacobianName), angleDeg, contrast, ...
        modeCount, clusterHC, clusterCosine, singularHC, singularCosine, ...
        norm(rightEnvelope), norm(leftEnvelope), ...
        norm(rightSingular), norm(leftSingular)};
end

result = cell2table(rows, 'VariableNames', { ...
    'jacobian', 'angleDeg', 'contrast', 'clusterModeCount', ...
    'normalizedEClusterHCnorm', 'normalizedEClusterCosine', ...
    'normalizedESingularHCnorm', 'normalizedESingularCosine', ...
    'rawRightEClusterNorm', 'rawLeftEClusterNorm', ...
    'rawRightESingularNorm', 'rawLeftESingularNorm'});
tag = sprintf('angle%s_contrast%d', ...
    strrep(sprintf('%.2f', angleDeg), '.', 'p'), contrast);
outputFile = fullfile(outputRoot, ...
    ['normalized_e_heatmap_left_right_' tag '.tsv']);
writetable(result, outputFile, ...
    'FileType', 'text', 'Delimiter', '\t');
mapFile = fullfile(outputRoot, ...
    ['normalized_e_heatmap_left_right_' tag '.mat']);
save(mapFile, 'conditionMaps', 'angleDeg', 'contrast', '-v7');
fprintf('Saved normalized E-map comparison to %s.\n', outputFile);
fprintf('Saved normalized E maps to %s.\n', mapFile);
end

function [rightEnvelope, leftEnvelope] = local_cluster_e_maps( ...
        matrix, matrixName, modeCount, n, eRows, iRows)
switch matrixName
    case 'J_full'
        active = matrix;
    case 'J_6'
        active = matrix(eRows, eRows);
    case 'J_I'
        active = matrix(iRows, iRows);
    otherwise
        error('NormalizedMap:Jacobian', 'Unsupported Jacobian %s.', matrixName);
end

rightActive = local_leading_vectors(active, modeCount);
leftActive = local_leading_vectors(active', modeCount);
rightVectors = complex(zeros(3*n, modeCount));
leftVectors = complex(zeros(3*n, modeCount));
switch matrixName
    case 'J_full'
        rightVectors = rightActive.Vectors;
        leftVectors = leftActive.Vectors;
    case 'J_6'
        for index = 1:modeCount
            lambda = rightActive.Values(index);
            rightVectors(:, index) = [rightActive.Vectors(:, index); ...
                matrix(iRows, eRows) * rightActive.Vectors(:, index) / lambda];
            leftVectors(:, index) = [leftActive.Vectors(:, index); ...
                zeros(n, 1)];
        end
    case 'J_I'
        for index = 1:modeCount
            lambda = rightActive.Values(index);
            rightVectors(:, index) = [ ...
                matrix(eRows, iRows) * rightActive.Vectors(:, index) / lambda; ...
                rightActive.Vectors(:, index)];
            leftVectors(:, index) = [zeros(2*n, 1); ...
                leftActive.Vectors(:, index)];
        end
end
rightBasis = local_basis(rightVectors);
leftBasis = local_basis(leftVectors);
rightEnvelope = local_e_envelope(rightBasis, n);
leftEnvelope = local_e_envelope(leftBasis, n);
end

function modes = local_leading_vectors(matrix, modeCount)
dimension = size(matrix, 1);
requested = min(dimension - 2, max(24, modeCount + 8));
options = struct('tol', 1e-10, 'maxit', 5000, ...
    'p', min(dimension, max(80, 2*requested + 16)), ...
    'disp', 0, 'isreal', true);
[vectors, diagonal, flag] = eigs( ...
    matrix, requested, 'largestreal', options);
if flag ~= 0
    error('NormalizedMap:Eigs', 'eigs returned flag %d.', flag);
end
values = diag(diagonal);
[~, order] = sortrows([real(values), imag(values)], [-1 -2]);
modes = struct('Vectors', vectors(:, order(1:modeCount)), ...
    'Values', values(order(1:modeCount)));
end

function basis = local_basis(vectors)
vectors = vectors ./ max(vecnorm(vectors), eps);
[basis, triangular] = qr(vectors, 0);
tolerance = max(size(triangular)) * eps(max(norm(triangular, 2), 1));
rankValue = sum(abs(diag(triangular)) > tolerance);
basis = basis(:, 1:rankValue);
end

function envelope = local_e_envelope(basis, n)
s = basis(1:n, :);
c = basis(n + (1:n), :);
e = 0.6923 * s + 0.3077 * c;
envelope = sqrt(sum(abs(e).^2, 2));
end

function [rightE, leftE] = local_singular_e_maps(matrix, n)
dimension = size(matrix, 1);
options = struct('tol', 1e-9, 'maxit', 3500, ...
    'p', min(dimension, 80), 'disp', 0);
[leftVector, ~, rightVector, flag] = ...
    svds(matrix, 1, 'largest', options);
if flag ~= 0
    error('NormalizedMap:Svds', 'svds returned flag %d.', flag);
end
phase = leftVector' * rightVector;
if abs(phase) > 0
    leftVector = leftVector * exp(1i * angle(phase));
end
rightE = real(0.6923 * rightVector(1:n) + ...
    0.3077 * rightVector(n + (1:n)));
leftE = real(0.6923 * leftVector(1:n) + ...
    0.3077 * leftVector(n + (1:n)));
end

function [hcNorm, cosine, rightMap, leftMap] = local_normalized_map_metrics( ...
        rightMap, leftMap, phaseInvariant)
rightMap = rightMap(:) / max(norm(rightMap(:)), eps);
leftMap = leftMap(:) / max(norm(leftMap(:)), eps);
if phaseInvariant && dot(rightMap, leftMap) < 0
    leftMap = -leftMap;
end
delta = rightMap - leftMap;
hcNorm = sqrt(0.8 * mean(abs(delta).^2));
cosine = real(rightMap' * leftMap);
if phaseInvariant
    cosine = abs(cosine);
end
end
