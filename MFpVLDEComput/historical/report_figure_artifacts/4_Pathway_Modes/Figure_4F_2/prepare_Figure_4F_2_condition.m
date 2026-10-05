function prepare_Figure_4F_2_condition(taskIndex, sourceRoot, modeSummaryFile, outputRoot)
% Prepare one condition for the order-independent Figure 4F.2 analysis.

% The finite-time analysis uses the Euclidean norm after observing the
% E-population map. Therefore W is identity in this 1600-dimensional E-map
% space (equivalently C_E' * C_E in the full S/C/I state space).

if nargin < 1 || isempty(taskIndex)
    taskIndex = str2double(getenv('SLURM_ARRAY_TASK_ID'));
end
if nargin < 2 || isempty(sourceRoot)
    sourceRoot = ['/scratch/xh2906/librarySCI_runs/' ...
        'h96_default_full_geometry_dedupe_array_localhelpers_20260607_084610/' ...
        'l6_mechanism_derivative_clamp_20260628_033413'];
end
if nargin < 3 || isempty(modeSummaryFile)
    modeSummaryFile = ['/scratch/xh2906/librarySCI_runs/' ...
        'h96_default_full_geometry_dedupe_array_localhelpers_20260607_084610/' ...
        'normalized_e_heatmap_left_right_20260807/code/' ...
        'left_right_mode_similarity_all_conditions.tsv'];
end
if nargin < 4 || isempty(outputRoot)
    outputRoot = fullfile(fileparts(mfilename('fullpath')), 'condition_data');
end
if ~isscalar(taskIndex) || taskIndex < 1 || taskIndex > 16 || ...
        taskIndex ~= round(taskIndex)
    error('Figure4F2:Task', 'Task index must be an integer from 1 to 16.');
end

angles = [0 7.5 15 22.5];
contrasts = [19 42 66 100];
[contrastIndex, angleIndex] = ind2sub([4 4], taskIndex);
angleDeg = angles(angleIndex);
contrast = contrasts(contrastIndex);

dynamicFile = local_geometry_file(sourceRoot, angleDeg, contrast, 0);
frozenFile = local_geometry_file(sourceRoot, angleDeg, contrast, 1);
dynamicData = load(dynamicFile, 'Section4', 'GeometryMetrics');
jFull = sparse(dynamicData.Section4.A);
taskTangents = double(dynamicData.Section4.TaskTangents);
featureNames = string(dynamicData.Section4.FeatureNames(:));
fixedPointE = double(dynamicData.GeometryMetrics.FixedPointMaps.E(:));
fullRightVectors = dynamicData.Section4.RightEigenVectors;
clear dynamicData

frozenData = load(frozenFile, 'Section4');
jFrozen = sparse(frozenData.Section4.A);
clear frozenData

dimension = size(jFull, 1);
n = dimension / 3;
if dimension ~= 4800 || n ~= round(n)
    error('Figure4F2:Dimension', 'Expected a 4800-dimensional Jacobian.');
end
eRows = 1:2*n;
iRows = 2*n + (1:n);
j6 = sparse(jFull - jFrozen);
jI = sparse(dimension, dimension);
jI(:, iRows) = jFull(:, iRows);
clear jFrozen

idxOrientation = find(featureNames == "orientation_deg", 1);
idxContrast = find(featureNames == "log_contrast", 1);
if isempty(idxOrientation) || isempty(idxContrast)
    error('Figure4F2:Features', 'Missing orientation or log-contrast tangent.');
end
qGain = fixedPointE;
qOrientation = local_e_component(taskTangents(:, idxOrientation), n);
qContrast = local_e_component(taskTangents(:, idxContrast), n);
Q = [qGain, qOrientation, qContrast];
gramMatrix = real(Q' * Q);
gramRcond = rcond(gramMatrix);
if gramRcond < 1e-12
    error('Figure4F2:VisualBasis', ...
        'Visual tangent Gram matrix is ill-conditioned (rcond %.3e).', gramRcond);
end
clear taskTangents

summary = readtable(modeSummaryFile, 'FileType', 'text', ...
    'Delimiter', '\t', 'TextType', 'string');
jacobianNames = ["J_full", "J_6", "J_I"];
modeCounts = zeros(1, 3);
for index = 1:3
    selected = summary.jacobian == jacobianNames(index) & ...
        abs(summary.angleDeg - angleDeg) < 1e-9 & ...
        summary.contrast == contrast;
    if sum(selected) ~= 1
        error('Figure4F2:ModeCount', 'Missing cluster count for %s.', ...
            jacobianNames(index));
    end
    modeCounts(index) = summary.clusterModeCount(selected);
end

preferredIndices = local_preferred_indices(angleDeg);
clusterVectors = cell(1, 3);
clusterVectors{1} = fullRightVectors(:, 1:modeCounts(1));
clusterVectors{2} = local_cluster_vectors( ...
    j6, 'J_6', modeCounts(2), n, eRows, iRows);
clusterVectors{3} = local_cluster_vectors( ...
    jI, 'J_I', modeCounts(3), n, eRows, iRows);
clear fullRightVectors

clusterFractions = zeros(3, 2);
singularFractions = zeros(3, 2);
clusterCoefficients = zeros(3, 3);
singularCoefficients = zeros(3, 3);
clusterShapley = zeros(3, 3);
singularShapley = zeros(3, 3);
matrices = {jFull, j6, jI};
singularValues = zeros(1, 3);
for index = 1:3
    signedClusterModes = local_signed_e_modes( ...
        clusterVectors{index}, n, preferredIndices);
    [modeFractions, modeCoefficients, modeShapley] = ...
        local_visual_fit(signedClusterModes, Q, gramMatrix);
    clusterFractions(index, :) = mean(modeFractions, 1);
    clusterCoefficients(index, :) = mean(modeCoefficients, 1);
    clusterShapley(index, :) = mean(modeShapley, 1);

    [leftVector, singularValues(index)] = ...
        local_top_singular_output(matrices{index});
    signedSingularMode = local_signed_e_modes( ...
        leftVector, n, preferredIndices);
    [modeFractions, modeCoefficients, modeShapley] = ...
        local_visual_fit(signedSingularMode, Q, gramMatrix);
    singularFractions(index, :) = modeFractions;
    singularCoefficients(index, :) = modeCoefficients;
    singularShapley(index, :) = modeShapley;
end

if ~exist(outputRoot, 'dir'); mkdir(outputRoot); end
outputFile = fullfile(outputRoot, sprintf( ...
    'Figure_4F_2_condition_angle%.2f_contrast%d.mat', angleDeg, contrast));
metricDescription = ['W = identity on the E observation map; full-state ' ...
    'equivalent W = C_E''*C_E with C_E = [0.6923I, 0.3077I, 0].'];
save(outputFile, 'taskIndex', 'angleDeg', 'contrast', ...
    'jacobianNames', 'modeCounts', 'clusterFractions', ...
    'singularFractions', 'clusterCoefficients', ...
    'singularCoefficients', 'clusterShapley', 'singularShapley', ...
    'singularValues', 'gramRcond', ...
    'metricDescription', '-v7');
fprintf('Saved Figure 4F.2 condition data to %s.\n', outputFile);
end

function fileName = local_geometry_file(sourceRoot, angleDeg, contrast, weight)
fileName = fullfile(sourceRoot, 'geometry_mat', sprintf( ...
    'geometry_sections4_5_h96baseline_L6eqW%dp00_contr%d_angle_%.2f.mat', ...
    weight, contrast, angleDeg));
end

function e = local_e_component(vector, n)
e = real(0.6923 * vector(1:n) + 0.3077 * vector(n + (1:n)));
end

function indices = local_preferred_indices(angleDeg)
if abs(angleDeg) < 1e-9
    pixel = [5 10];
elseif abs(angleDeg - 7.5) < 1e-9
    pixel = [4 10];
elseif abs(angleDeg - 15) < 1e-9
    pixel = [3 10];
elseif abs(angleDeg - 22.5) < 1e-9
    pixel = [1 10];
else
    error('Figure4F2:Angle', 'Unsupported angle %.2f.', angleDeg);
end
rows = pixel(1):10:(pixel(1) + 30);
columns = pixel(2):10:(pixel(2) + 30);
[columnGrid, rowGrid] = meshgrid(columns, rows);
indices = sub2ind([40 40], rowGrid(:), columnGrid(:));
end

function vectors = local_cluster_vectors(matrix, matrixName, modeCount, n, eRows, iRows)
if strcmp(matrixName, 'J_6')
    active = matrix(eRows, eRows);
else
    active = matrix(iRows, iRows);
end
options = struct('tol', 1e-10, 'maxit', 5000, ...
    'p', min(size(active, 1), max(80, 2*modeCount + 16)), ...
    'disp', 0, 'isreal', true);
[activeVectors, diagonal, flag] = eigs( ...
    active, modeCount, 'largestreal', options);
if flag ~= 0
    error('Figure4F2:Eigs', 'eigs returned flag %d for %s.', flag, matrixName);
end
values = diag(diagonal);
[~, order] = sortrows([real(values), imag(values)], [-1 -2]);
activeVectors = activeVectors(:, order);
values = values(order);
vectors = complex(zeros(3*n, modeCount));
for index = 1:modeCount
    if strcmp(matrixName, 'J_6')
        vectors(:, index) = [activeVectors(:, index); ...
            matrix(iRows, eRows) * activeVectors(:, index) / values(index)];
    else
        vectors(:, index) = [ ...
            matrix(eRows, iRows) * activeVectors(:, index) / values(index); ...
            activeVectors(:, index)];
    end
end
end

function signedModes = local_signed_e_modes(vectors, n, preferredIndices)
signedModes = zeros(n, size(vectors, 2));
for index = 1:size(vectors, 2)
    eComplex = 0.6923 * vectors(1:n, index) + ...
        0.3077 * vectors(n + (1:n), index);
    anchor = sum(eComplex(preferredIndices));
    if abs(anchor) < 1e-12
        [~, anchorIndex] = max(abs(eComplex));
        anchor = eComplex(anchorIndex);
    end
    e = real(eComplex * exp(-1i * angle(anchor)));
    if sum(e(preferredIndices)) < 0; e = -e; end
    signedModes(:, index) = e / max(norm(e), eps);
end
end

function [leftVector, singularValue] = local_top_singular_output(matrix)
options = struct('tol', 1e-9, 'maxit', 5000, ...
    'p', min(size(matrix, 1), 80), 'disp', 0);
[leftVector, diagonal, ~, flag] = svds(matrix, 1, 'largest', options);
if flag ~= 0; error('Figure4F2:Svds', 'svds returned flag %d.', flag); end
singularValue = diagonal(1, 1);
end

function [fractions, coefficients, shapley] = local_visual_fit(modes, Q, gramMatrix)
modeCount = size(modes, 2);
fractions = zeros(modeCount, 2);
coefficients = zeros(modeCount, 3);
shapley = zeros(modeCount, 3);
for index = 1:modeCount
    mode = modes(:, index);
    mode = mode / max(norm(mode), eps);
    coefficient = gramMatrix \ real(Q' * mode);
    projected = Q * coefficient;
    visualFraction = real(projected' * projected) / real(mode' * mode);
    visualFraction = min(max(visualFraction, 0), 1);
    fractions(index, :) = [visualFraction, 1 - visualFraction];
    coefficients(index, :) = coefficient(:).';
    shapley(index, :) = local_shapley_partition(mode, Q, visualFraction);
end
end

function phi = local_shapley_partition(mode, Q, fullVisualFraction)
% Order-independent variance partition across gain, orientation, contrast.
featureCount = size(Q, 2);
subsetR2 = zeros(1, 2^featureCount);
for mask = 1:(2^featureCount - 1)
    selected = logical(bitget(mask, 1:featureCount));
    subsetQ = Q(:, selected);
    subsetGram = real(subsetQ' * subsetQ);
    subsetCoefficient = subsetGram \ real(subsetQ' * mode);
    subsetProjection = subsetQ * subsetCoefficient;
    subsetR2(mask + 1) = real(subsetProjection' * subsetProjection) / ...
        real(mode' * mode);
end

phi = zeros(1, featureCount);
for feature = 1:featureCount
    for mask = 0:(2^featureCount - 1)
        if bitget(mask, feature); continue; end
        subsetSize = sum(bitget(mask, 1:featureCount));
        weight = factorial(subsetSize) * ...
            factorial(featureCount - subsetSize - 1) / factorial(featureCount);
        withFeature = bitset(mask, feature, 1);
        phi(feature) = phi(feature) + weight * ...
            (subsetR2(withFeature + 1) - subsetR2(mask + 1));
    end
end
if any(phi < -1e-10)
    error('Figure4F2:Shapley', ...
        'Negative Shapley contribution exceeds numerical tolerance.');
end
phi = max(phi, 0);
if fullVisualFraction > eps
    phi = phi * (fullVisualFraction / sum(phi));
else
    phi(:) = 0;
end
if abs(sum(phi) - fullVisualFraction) > 1e-10
    error('Figure4F2:ShapleySum', 'Shapley values do not sum to visual R^2.');
end
end
