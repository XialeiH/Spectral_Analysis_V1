function conditionData = run_full_eigenspectrum_maps_condition_task( ...
        taskIndex, sourceRoot, eigenvalueCacheRoot, outputRoot)
% Extract compact map and spectrum data for one of the 16 requested conditions.

% Matrix definitions are kept identical to the corrected Jacobian analyses:
%   J_full = baseline dynamic-L6 Jacobian
%   J_6    = [D_*l D_lE,S^dyn, D_*l D_lE,C^dyn, 0]
%          = J_full - frozen-L6 Jacobian at the same fixed point
%   J_I    = inhibitory-source columns of J_full

% The top-real cluster is represented by the coordinate envelope of its
% orthogonal projector. This is invariant to eigenvector mixing inside a
% repeated or tightly clustered eigenspace.

if nargin < 1 || isempty(taskIndex)
    taskIndex = str2double(getenv('FULLMAP_TASK_INDEX'));
    if ~isfinite(taskIndex)
        taskIndex = str2double(getenv('SLURM_ARRAY_TASK_ID'));
    end
end
if nargin < 2 || isempty(sourceRoot)
    sourceRoot = getenv('FULLMAP_SOURCE_ROOT');
end
if nargin < 3 || isempty(eigenvalueCacheRoot)
    eigenvalueCacheRoot = getenv('FULLMAP_EIGEN_CACHE_ROOT');
end
if nargin < 4 || isempty(outputRoot)
    outputRoot = getenv('FULLMAP_OUTPUT_ROOT');
end

angles = [0 7.5 15 22.5];
contrasts = [19 42 66 100];
conditionCount = numel(angles) * numel(contrasts);
if ~isscalar(taskIndex) || ~isfinite(taskIndex) || taskIndex ~= round(taskIndex) || ...
        taskIndex < 1 || taskIndex > conditionCount
    error('FullMaps:TaskIndex', 'Task index must be an integer from 1 to 16.');
end
if ~isfolder(sourceRoot)
    error('FullMaps:SourceRoot', 'Missing source root: %s', sourceRoot);
end
if ~isfolder(eigenvalueCacheRoot)
    error('FullMaps:EigenCache', ...
        'Missing compact eigenvalue cache: %s', eigenvalueCacheRoot);
end
if isempty(outputRoot)
    error('FullMaps:OutputRoot', 'An output root is required.');
end
if ~exist(outputRoot, 'dir'); mkdir(outputRoot); end
dataRoot = fullfile(outputRoot, 'data');
if ~exist(dataRoot, 'dir'); mkdir(dataRoot); end

[contrastIndex, angleIndex] = ind2sub( ...
    [numel(contrasts), numel(angles)], taskIndex);
angleValue = angles(angleIndex);
contrast = contrasts(contrastIndex);

file0 = fullfile(sourceRoot, 'geometry_mat', sprintf( ...
    'geometry_sections4_5_h96baseline_L6eqW0p00_contr%d_angle_%.2f.mat', ...
    contrast, angleValue));
file1 = fullfile(sourceRoot, 'geometry_mat', sprintf( ...
    'geometry_sections4_5_h96baseline_L6eqW1p00_contr%d_angle_%.2f.mat', ...
    contrast, angleValue));
local_require_file(file0, 'dynamic-L6 geometry');
local_require_file(file1, 'frozen-L6 geometry');

loaded0 = load(file0, 'Section4', 'GeometryMetrics');
loaded1 = load(file1, 'Section4');
if ~isfield(loaded0, 'Section4') || ~isfield(loaded0.Section4, 'A') || ...
        ~isfield(loaded1, 'Section4') || ~isfield(loaded1.Section4, 'A')
    error('FullMaps:MissingJacobian', 'Section4.A is missing from a source file.');
end
if ~isfield(loaded0, 'GeometryMetrics') || ...
        ~isfield(loaded0.GeometryMetrics, 'FixedPointMaps')
    error('FullMaps:MissingFixedPoint', ...
        'GeometryMetrics.FixedPointMaps is missing from %s.', file0);
end

jFull = sparse(loaded0.Section4.A);
jNoL6 = sparse(loaded1.Section4.A);
dimension = size(jFull, 1);
if size(jFull, 2) ~= dimension || dimension ~= 4800 || ...
        ~isequal(size(jFull), size(jNoL6))
    error('FullMaps:JacobianSize', ...
        'Expected two compatible 4800-by-4800 Jacobians.');
end
n = dimension / 3;
side = round(sqrt(n));
if side^2 ~= n
    error('FullMaps:MapSize', 'Population size %d is not a square map.', n);
end

cacheFile = fullfile(eigenvalueCacheRoot, sprintf( ...
    'near_zero_eigenvalues_angle%s_contrast%d.mat', ...
    local_angle_tag(angleValue), contrast));
local_require_file(cacheFile, 'compact eigenvalue cache');
cache = load(cacheFile, 'eigenvalueData');
requiredSpectra = {'J_baseline', 'J_6', 'J_I'};
for fieldIndex = 1:numel(requiredSpectra)
    if ~isfield(cache.eigenvalueData, requiredSpectra{fieldIndex})
        error('FullMaps:EigenCacheField', ...
            'Missing eigenvalueData.%s in %s.', ...
            requiredSpectra{fieldIndex}, cacheFile);
    end
end

fixedPointMaps = local_fixed_point_maps( ...
    loaded0.GeometryMetrics.FixedPointMaps, side);
clear loaded0 loaded1

j6 = sparse(jFull - jNoL6);
clear jNoL6
iColumns = 2 * n + (1:n);
jI = sparse(dimension, dimension);
jI(:, iColumns) = jFull(:, iColumns);

decompositionAudit = struct();
decompositionAudit.J6ISourceFraction = ...
    norm(j6(:, iColumns), 'fro') / max(norm(j6, 'fro'), eps);
decompositionAudit.JINonISourceFraction = ...
    norm(jI(:, 1:2*n), 'fro') / max(norm(jI, 'fro'), eps);
if decompositionAudit.J6ISourceFraction > 1e-12 || ...
        decompositionAudit.JINonISourceFraction > 1e-12
    error('FullMaps:Decomposition', ...
        'The source-column decomposition failed its structural audit.');
end

wC = 0.3077;
clusterTolerance = 1e-3;
analyses = struct();
analyses.J_full = local_matrix_analysis(jFull, 'J_full', ...
    cache.eigenvalueData.J_baseline, n, side, wC, clusterTolerance);
analyses.J_6 = local_matrix_analysis(j6, 'J_6', ...
    cache.eigenvalueData.J_6, n, side, wC, clusterTolerance);
analyses.J_I = local_matrix_analysis(jI, 'J_I', ...
    cache.eigenvalueData.J_I, n, side, wC, clusterTolerance);

conditionData = struct();
conditionData.AngleDeg = angleValue;
conditionData.Contrast = contrast;
conditionData.TaskIndex = taskIndex;
conditionData.PopulationSize = n;
conditionData.MapSide = side;
conditionData.CWeight = wC;
conditionData.ClusterTolerance = clusterTolerance;
conditionData.FixedPointMaps = fixedPointMaps;
conditionData.Jacobians = analyses;
conditionData.DecompositionAudit = decompositionAudit;
conditionData.SourceFiles = struct('DynamicL6', file0, ...
    'FrozenL6', file1, 'EigenvalueCache', cacheFile);

tag = sprintf('angle%s_contrast%d', local_angle_tag(angleValue), contrast);
outputFile = fullfile(dataRoot, ['full_eigenspectrum_maps_', tag, '.mat']);
save(outputFile, 'conditionData', '-v7');

summaryRows = cell(3, 10);
matrixNames = {'J_full', 'J_6', 'J_I'};
for matrixIndex = 1:numel(matrixNames)
    name = matrixNames{matrixIndex};
    item = analyses.(name);
    summaryRows(matrixIndex, :) = {string(name), angleValue, contrast, ...
        max(real(item.Eigenvalues)), min(real(item.Eigenvalues)), ...
        max(abs(item.Eigenvalues)), item.TopEigenvalue, ...
        item.TopClusterCount, item.TopClusterRank, item.LeadingSingularValue};
end
summary = cell2table(summaryRows, 'VariableNames', { ...
    'jacobian', 'angleDeg', 'contrast', 'maximumRealEigenvalue', ...
    'minimumRealEigenvalue', 'spectralRadius', 'topEigenvalue', ...
    'topClusterCount', 'topClusterRank', 'leadingSingularValue'});
writetable(summary, fullfile(dataRoot, ...
    ['full_eigenspectrum_maps_summary_', tag, '.tsv']), ...
    'FileType', 'text', 'Delimiter', '\t');

fprintf(['Full eigenspectrum map extraction complete: angle %.2f deg, ' ...
    'contrast %d -> %s\n'], angleValue, contrast, outputFile);
end

function analysis = local_matrix_analysis(matrix, matrixName, eigenvalues, ...
        n, side, wC, clusterTolerance)
eigenvalues = eigenvalues(:);
if numel(eigenvalues) ~= 3 * n
    error('FullMaps:EigenvalueCount', ...
        '%s has %d cached eigenvalues; expected %d.', ...
        matrixName, numel(eigenvalues), 3*n);
end

[topVector, topEigenvalue, clusterBasis, clusterEigenvalues, ...
    topResidual] = local_top_real_modes(matrix, matrixName, eigenvalues, ...
    n, clusterTolerance);
[singularOutput, singularValue, singularResidual] = ...
    local_leading_singular_output(matrix);

analysis = struct();
analysis.Eigenvalues = eigenvalues;
analysis.TopEigenvalue = topEigenvalue;
analysis.TopEigenvectorMaps = local_vector_maps(topVector, n, side, wC);
analysis.TopClusterEnvelopeMaps = ...
    local_cluster_envelope_maps(clusterBasis, n, side, wC);
analysis.TopClusterEigenvalues = clusterEigenvalues;
analysis.TopClusterCount = numel(clusterEigenvalues);
analysis.TopClusterRank = size(clusterBasis, 2);
analysis.TopEigenResidual = topResidual;
analysis.LeadingSingularOutputMaps = ...
    local_vector_maps(singularOutput, n, side, wC);
analysis.LeadingSingularValue = singularValue;
analysis.LeadingSingularResidual = singularResidual;
end

function [topVector, topEigenvalue, clusterBasis, clusterEigenvalues, ...
        residual] = local_top_real_modes(matrix, matrixName, eigenvalues, ...
        n, clusterTolerance)
dimension = size(matrix, 1);
eRows = 1:2*n;
iRows = 2*n + (1:n);
switch matrixName
    case 'J_full'
        activeMatrix = matrix;
        liftMode = 'identity';
    case 'J_6'
        activeMatrix = matrix(eRows, eRows);
        liftMode = 'J6';
    case 'J_I'
        activeMatrix = matrix(iRows, iRows);
        liftMode = 'JI';
    otherwise
        error('FullMaps:MatrixName', 'Unsupported matrix %s.', matrixName);
end

topReal = max(real(eigenvalues));
selectionTolerance = clusterTolerance + ...
    100 * eps(max(abs(topReal), 1));
expectedCount = sum(real(eigenvalues) >= topReal - selectionTolerance);
activeDimension = size(activeMatrix, 1);
requestedCount = min(activeDimension - 2, max(16, expectedCount + 8));
maximumCount = min(activeDimension - 2, max(192, expectedCount + 32));

while true
    options = struct('tol', 1e-10, 'maxit', 5000, ...
        'p', min(activeDimension, max(80, 2*requestedCount + 16)), ...
        'disp', 0, 'isreal', true);
    [activeVectors, diagonal, flag] = eigs( ...
        activeMatrix, requestedCount, 'largestreal', options);
    if flag ~= 0
        error('FullMaps:EigsConvergence', ...
            'eigs returned flag %d for %s.', flag, matrixName);
    end
    activeValues = diag(diagonal);
    [~, order] = sortrows([real(activeValues), imag(activeValues)], [-1 -2]);
    activeValues = activeValues(order);
    activeVectors = activeVectors(:, order);
    selected = find(real(activeValues) >= topReal - selectionTolerance);
    if numel(selected) >= expectedCount || requestedCount >= maximumCount
        break
    end
    requestedCount = min(maximumCount, max(2*requestedCount, expectedCount + 8));
end
if numel(selected) < expectedCount
    error('FullMaps:ClusterIncomplete', ...
        '%s recovered %d of %d top-cluster modes.', ...
        matrixName, numel(selected), expectedCount);
end
selected = selected(1:expectedCount);
clusterEigenvalues = activeValues(selected);

lifted = complex(zeros(dimension, expectedCount));
for modeIndex = 1:expectedCount
    vector = activeVectors(:, selected(modeIndex));
    lambda = clusterEigenvalues(modeIndex);
    switch liftMode
        case 'identity'
            fullVector = vector;
        case 'J6'
            if abs(lambda) <= 1e-12
                error('FullMaps:J6ZeroMode', ...
                    'The selected J_6 top cluster unexpectedly contains zero.');
            end
            fullVector = [vector; matrix(iRows, eRows) * vector / lambda];
        case 'JI'
            if abs(lambda) <= 1e-12
                error('FullMaps:JIZeroMode', ...
                    'The selected J_I top cluster unexpectedly contains zero.');
            end
            fullVector = [matrix(eRows, iRows) * vector / lambda; vector];
    end
    lifted(:, modeIndex) = fullVector / max(norm(fullVector), eps);
end

topVector = local_canonical_phase(lifted(:, 1));
topEigenvalue = clusterEigenvalues(1);
residual = norm(matrix * topVector - topEigenvalue * topVector) / ...
    max((norm(matrix, 'fro') + abs(topEigenvalue)) * norm(topVector), eps);

[basis, triangularFactor] = qr(lifted, 0);
rankTolerance = max(size(triangularFactor)) * ...
    eps(max(norm(triangularFactor, 2), 1));
clusterRank = sum(abs(diag(triangularFactor)) > rankTolerance);
clusterBasis = basis(:, 1:clusterRank);
end

function [outputVector, singularValue, residual] = ...
        local_leading_singular_output(matrix)
dimension = size(matrix, 1);
options = struct('tol', 1e-9, 'maxit', 3500, ...
    'p', min(dimension, 80), 'disp', 0);
[outputVector, singularMatrix, inputVector, flag] = ...
    svds(matrix, 1, 'largest', options);
if flag ~= 0
    error('FullMaps:SvdsConvergence', 'svds returned flag %d.', flag);
end
singularValue = singularMatrix(1, 1);
outputVector = local_canonical_phase(outputVector);
phase = outputVector' * (matrix * inputVector);
if real(phase) < 0
    inputVector = -inputVector;
end
residual = norm(matrix * inputVector - singularValue * outputVector) / ...
    max(norm(matrix, 'fro'), eps);
end

function maps = local_fixed_point_maps(sourceMaps, side)
required = {'S', 'C', 'I', 'E'};
maps = struct();
for index = 1:numel(required)
    name = required{index};
    if ~isfield(sourceMaps, name) || numel(sourceMaps.(name)) ~= side^2
        error('FullMaps:FixedPointMap', ...
            'FixedPointMaps.%s is missing or has the wrong size.', name);
    end
    maps.(name) = reshape(real(sourceMaps.(name)), side, side);
end
end

function maps = local_vector_maps(vector, n, side, wC)
wS = 1 - wC;
s = vector(1:n);
c = vector(n + (1:n));
i = vector(2*n + (1:n));
e = wS*s + wC*c;
maps = struct('S', reshape(real(s), side, side), ...
    'C', reshape(real(c), side, side), ...
    'I', reshape(real(i), side, side), ...
    'E', reshape(real(e), side, side));
end

function maps = local_cluster_envelope_maps(basis, n, side, wC)
wS = 1 - wC;
s = basis(1:n, :);
c = basis(n + (1:n), :);
i = basis(2*n + (1:n), :);
e = wS*s + wC*c;
maps = struct('S', reshape(sqrt(sum(abs(s).^2, 2)), side, side), ...
    'C', reshape(sqrt(sum(abs(c).^2, 2)), side, side), ...
    'I', reshape(sqrt(sum(abs(i).^2, 2)), side, side), ...
    'E', reshape(sqrt(sum(abs(e).^2, 2)), side, side));
end

function vector = local_canonical_phase(vector)
[~, pivot] = max(abs(vector));
if abs(vector(pivot)) > 0
    vector = vector * exp(-1i * angle(vector(pivot)));
end
if real(vector(pivot)) < 0
    vector = -vector;
end
vector = vector / max(norm(vector), eps);
end

function local_require_file(pathValue, description)
if ~isfile(pathValue)
    error('FullMaps:MissingFile', 'Missing %s: %s', description, pathValue);
end
end

function tag = local_angle_tag(value)
tag = strrep(sprintf('%.2f', value), '.', 'p');
end
