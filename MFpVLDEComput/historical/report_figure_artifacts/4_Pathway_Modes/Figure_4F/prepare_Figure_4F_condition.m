function prepare_Figure_4F_condition(taskIndex, sourceRoot, modeSummaryFile, outputRoot)
% Prepare one condition of signed visual-geometry projections for Figure 4F.

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
    error('Figure4F:Task', 'Task index must be an integer from 1 to 16.');
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
taskTangents = dynamicData.Section4.TaskTangents;
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
    error('Figure4F:Dimension', 'Expected a 4800-dimensional Jacobian.');
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
    error('Figure4F:Features', 'Missing task tangents.');
end
qOrientation = local_e_component(taskTangents(:, idxOrientation), n);
qContrast = local_e_component(taskTangents(:, idxContrast), n);
visualBasis = local_ordered_basis([fixedPointE, qOrientation, qContrast]);
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
        error('Figure4F:ModeCount', 'Missing cluster count for %s.', ...
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

clusterProjections = zeros(3, 4);
singularProjections = zeros(3, 4);
matrices = {jFull, j6, jI};
singularValues = zeros(1, 3);
for index = 1:3
    signedClusterModes = local_signed_e_modes( ...
        clusterVectors{index}, n, preferredIndices);
    clusterProjections(index, :) = mean( ...
        local_project_modes(signedClusterModes, visualBasis), 1);
    [leftVector, singularValues(index)] = ...
        local_top_singular_output(matrices{index});
    signedSingularMode = local_signed_e_modes( ...
        leftVector, n, preferredIndices);
    singularProjections(index, :) = ...
        local_project_modes(signedSingularMode, visualBasis);
end

if ~exist(outputRoot, 'dir'); mkdir(outputRoot); end
outputFile = fullfile(outputRoot, sprintf( ...
    'Figure_4F_condition_angle%.2f_contrast%d.mat', angleDeg, contrast));
save(outputFile, 'taskIndex', 'angleDeg', 'contrast', ...
    'jacobianNames', 'modeCounts', 'clusterProjections', ...
    'singularProjections', 'singularValues', '-v7');
fprintf('Saved Figure 4F condition data to %s.\n', outputFile);
end

function fileName = local_geometry_file(sourceRoot, angleDeg, contrast, weight)
fileName = fullfile(sourceRoot, 'geometry_mat', sprintf( ...
    'geometry_sections4_5_h96baseline_L6eqW%dp00_contr%d_angle_%.2f.mat', ...
    weight, contrast, angleDeg));
end

function e = local_e_component(vector, n)
e = real(0.6923 * vector(1:n) + 0.3077 * vector(n + (1:n)));
end

function basis = local_ordered_basis(directions)
basis = zeros(size(directions));
for column = 1:size(directions, 2)
    vector = directions(:, column);
    for previous = 1:(column - 1)
        vector = vector - basis(:, previous) * (basis(:, previous)' * vector);
    end
    vectorNorm = norm(vector);
    if vectorNorm < 1e-10 * max(norm(directions(:, column)), 1)
        error('Figure4F:VisualBasis', 'Dependent visual direction %d.', column);
    end
    basis(:, column) = vector / vectorNorm;
end
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
    error('Figure4F:Angle', 'Unsupported angle %.2f.', angleDeg);
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
    error('Figure4F:Eigs', 'eigs returned flag %d for %s.', flag, matrixName);
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
if flag ~= 0; error('Figure4F:Svds', 'svds returned flag %d.', flag); end
singularValue = diagonal(1, 1);
end

function projections = local_project_modes(modes, basis)
projections = zeros(size(modes, 2), size(basis, 2) + 1);
for index = 1:size(modes, 2)
    coefficients = basis' * modes(:, index);
    residual = modes(:, index) - basis * coefficients;
    projections(index, :) = [coefficients(:).', norm(residual)];
end
end
