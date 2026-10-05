function prepare_Figure_4F_data(sourceRoot, outputFile)
% Prepare signed visual-geometry projections for Figure 4F.

if nargin < 1 || isempty(sourceRoot)
    sourceRoot = ['/scratch/xh2906/librarySCI_runs/' ...
        'h96_default_full_geometry_dedupe_array_localhelpers_20260607_084610/' ...
        'l6_mechanism_derivative_clamp_20260628_033413'];
end
if nargin < 2 || isempty(outputFile)
    outputFile = fullfile(fileparts(mfilename('fullpath')), ...
        'Figure_4F_projection_data.mat');
end

angleDeg = 0;
contrast = 100;
dynamicFile = local_geometry_file(sourceRoot, angleDeg, contrast, 0);
frozenFile = local_geometry_file(sourceRoot, angleDeg, contrast, 1);

dynamicData = load(dynamicFile, 'Section4', 'GeometryMetrics');
jFull = sparse(dynamicData.Section4.A);
taskTangents = dynamicData.Section4.TaskTangents;
featureNames = string(dynamicData.Section4.FeatureNames(:));
fixedPointE = double(dynamicData.GeometryMetrics.FixedPointMaps.E(:));
topFullVector = dynamicData.Section4.RightEigenVectors(:, 1);
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
    error('Figure4F:Features', ...
        'Expected orientation_deg and log_contrast task tangents.');
end
qOrientation = local_e_component(taskTangents(:, idxOrientation), n);
qContrast = local_e_component(taskTangents(:, idxContrast), n);
clear taskTangents

rawDirections = [fixedPointE, qOrientation, qContrast];
[visualBasis, directionNorms] = local_ordered_basis(rawDirections);
componentNames = ["Gain", "Orientation", "Contrast", ...
    "Orthogonal residual"];
jacobianNames = ["J_full", "J_6", "J_I"];

preferredIndices = local_preferred_indices();
eigenModes = zeros(n, 3);
singularModes = zeros(n, 3);
leadingEigenvalues = complex(zeros(1, 3));
singularValues = zeros(1, 3);

eigenModes(:, 1) = local_signed_e_mode(topFullVector, n, preferredIndices);
leadingEigenvalues(1) = local_rayleigh(jFull, topFullVector);

[vector6, leadingEigenvalues(2)] = local_leading_right_mode( ...
    j6, 'J_6', n, eRows, iRows);
eigenModes(:, 2) = local_signed_e_mode(vector6, n, preferredIndices);
[vectorI, leadingEigenvalues(3)] = local_leading_right_mode( ...
    jI, 'J_I', n, eRows, iRows);
eigenModes(:, 3) = local_signed_e_mode(vectorI, n, preferredIndices);

matrices = {jFull, j6, jI};
for index = 1:3
    [leftVector, singularValues(index)] = local_top_singular_output(matrices{index});
    singularModes(:, index) = local_signed_e_mode( ...
        leftVector, n, preferredIndices);
end

eigenProjections = local_project_modes(eigenModes, visualBasis);
singularProjections = local_project_modes(singularModes, visualBasis);

projectionTable = table();
for modeIndex = 1:2
    if modeIndex == 1
        modeName = "Leading eigenmode";
        values = eigenProjections;
    else
        modeName = "Finite-time singular response";
        values = singularProjections;
    end
    for jacobianIndex = 1:3
        next = table(repmat(modeName, 4, 1), ...
            repmat(jacobianNames(jacobianIndex), 4, 1), ...
            componentNames(:), values(jacobianIndex, :).', ...
            'VariableNames', {'Mode', 'Jacobian', 'Component', 'Projection'});
        projectionTable = [projectionTable; next]; %#ok<AGROW>
    end
end

outputFolder = fileparts(outputFile);
if ~exist(outputFolder, 'dir'); mkdir(outputFolder); end
save(outputFile, 'angleDeg', 'contrast', ...
    'componentNames', 'jacobianNames', 'visualBasis', 'directionNorms', ...
    'eigenModes', 'singularModes', 'eigenProjections', ...
    'singularProjections', 'leadingEigenvalues', 'singularValues', ...
    'projectionTable', '-v7');
writetable(projectionTable, strrep(outputFile, '.mat', '.tsv'), ...
    'FileType', 'text', 'Delimiter', '\t');
fprintf('Saved Figure 4F projection data to %s.\n', outputFile);
disp(projectionTable);
end

function fileName = local_geometry_file(sourceRoot, angleDeg, contrast, weight)
fileName = fullfile(sourceRoot, 'geometry_mat', sprintf( ...
    'geometry_sections4_5_h96baseline_L6eqW%d%s_contr%d_angle_%.2f.mat', ...
    floor(weight), local_weight_suffix(weight), contrast, angleDeg));
end

function suffix = local_weight_suffix(weight)
if weight == 0
    suffix = 'p00';
elseif weight == 1
    suffix = 'p00';
else
    error('Figure4F:Weight', 'Only weights 0 and 1 are supported.');
end
end

function e = local_e_component(vector, n)
e = real(0.6923 * vector(1:n) + 0.3077 * vector(n + (1:n)));
end

function [basis, rawNorms] = local_ordered_basis(directions)
rawNorms = vecnorm(directions, 2, 1);
basis = zeros(size(directions));
for column = 1:size(directions, 2)
    vector = directions(:, column);
    for previous = 1:(column - 1)
        vector = vector - basis(:, previous) * (basis(:, previous)' * vector);
    end
    vectorNorm = norm(vector);
    if vectorNorm < 1e-10 * max(rawNorms(column), 1)
        error('Figure4F:VisualBasis', ...
            'Visual direction %d is linearly dependent.', column);
    end
    basis(:, column) = vector / vectorNorm;
end
end

function indices = local_preferred_indices()
% Angle-0 preferred pixel is (5,10) in each 10-by-10 HC block.
rows = 5:10:35;
columns = 10:10:40;
[columnGrid, rowGrid] = meshgrid(columns, rows);
indices = sub2ind([40 40], rowGrid(:), columnGrid(:));
end

function signedMode = local_signed_e_mode(vector, n, preferredIndices)
eComplex = 0.6923 * vector(1:n) + 0.3077 * vector(n + (1:n));
anchor = sum(eComplex(preferredIndices));
if abs(anchor) < 1e-12
    [~, anchorIndex] = max(abs(eComplex));
    anchor = eComplex(anchorIndex);
end
eComplex = eComplex * exp(-1i * angle(anchor));
signedMode = real(eComplex);
if sum(signedMode(preferredIndices)) < 0
    signedMode = -signedMode;
end
signedMode = signedMode / max(norm(signedMode), eps);
end

function lambda = local_rayleigh(matrix, vector)
lambda = (vector' * (matrix * vector)) / (vector' * vector);
end

function [vector, lambda] = local_leading_right_mode(matrix, matrixName, n, eRows, iRows)
switch matrixName
    case 'J_6'
        active = matrix(eRows, eRows);
    case 'J_I'
        active = matrix(iRows, iRows);
    otherwise
        error('Figure4F:Jacobian', 'Unsupported Jacobian %s.', matrixName);
end
options = struct('tol', 1e-10, 'maxit', 5000, ...
    'p', min(size(active, 1), 80), 'disp', 0, 'isreal', true);
[activeVector, diagonal, flag] = eigs(active, 1, 'largestreal', options);
if flag ~= 0
    error('Figure4F:Eigs', 'eigs returned flag %d for %s.', flag, matrixName);
end
lambda = diagonal(1, 1);
switch matrixName
    case 'J_6'
        vector = [activeVector; matrix(iRows, eRows) * activeVector / lambda];
    case 'J_I'
        vector = [matrix(eRows, iRows) * activeVector / lambda; activeVector];
end
end

function [leftVector, singularValue] = local_top_singular_output(matrix)
options = struct('tol', 1e-9, 'maxit', 5000, ...
    'p', min(size(matrix, 1), 80), 'disp', 0);
[leftVector, diagonal, ~, flag] = svds(matrix, 1, 'largest', options);
if flag ~= 0
    error('Figure4F:Svds', 'svds returned flag %d.', flag);
end
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
