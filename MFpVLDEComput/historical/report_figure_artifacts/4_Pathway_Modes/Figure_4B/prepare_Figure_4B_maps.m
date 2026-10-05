function prepare_Figure_4B_maps(sourceRoot, modeSummaryFile, outputFile)
% Compute normalized E/I right and left leading-cluster maps for Figure 4B.

if nargin < 1 || isempty(sourceRoot)
    sourceRoot = ['/scratch/xh2906/librarySCI_runs/' ...
        'h96_default_full_geometry_dedupe_array_localhelpers_20260607_084610/' ...
        'l6_mechanism_derivative_clamp_20260628_033413'];
end
if nargin < 2 || isempty(modeSummaryFile)
    modeSummaryFile = ['/scratch/xh2906/librarySCI_runs/' ...
        'h96_default_full_geometry_dedupe_array_localhelpers_20260607_084610/' ...
        'normalized_e_heatmap_left_right_20260807/code/' ...
        'left_right_mode_similarity_all_conditions.tsv'];
end
if nargin < 3 || isempty(outputFile)
    outputFile = fullfile(fileparts(mfilename('fullpath')), ...
        'Figure_4B_cluster_maps.mat');
end

angleDeg = 0;
contrast = 100;
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
    error('Figure4B:Dimension', 'Expected a 4800-dimensional Jacobian.');
end
eRows = 1:2*n;
iRows = 2*n + (1:n);
j6 = sparse(jFull - jFrozen);
jI = sparse(dimension, dimension);
jI(:, iRows) = jFull(:, iRows);

modeSummary = readtable(modeSummaryFile, ...
    'FileType', 'text', 'Delimiter', '\t', 'TextType', 'string');
jacobianNames = {'J_full', 'J_6', 'J_I'};
matrices = {jFull, j6, jI};
clusterMaps = struct();
modeCounts = zeros(1, 3);
for jacobianIndex = 1:3
    jacobianName = jacobianNames{jacobianIndex};
    selected = modeSummary.jacobian == string(jacobianName) & ...
        abs(modeSummary.angleDeg - angleDeg) < 1e-9 & ...
        modeSummary.contrast == contrast;
    if sum(selected) ~= 1
        error('Figure4B:ModeCount', ...
            'Missing unique mode count for %s.', jacobianName);
    end
    modeCounts(jacobianIndex) = modeSummary.clusterModeCount(selected);
    [rightBasis, leftBasis] = local_cluster_bases( ...
        matrices{jacobianIndex}, jacobianName, ...
        modeCounts(jacobianIndex), n, eRows, iRows);
    clusterMaps.(jacobianName).ERight = reshape( ...
        local_normalize(local_e_envelope(rightBasis, n)), 40, 40);
    clusterMaps.(jacobianName).ELeft = reshape( ...
        local_normalize(local_e_envelope(leftBasis, n)), 40, 40);
    clusterMaps.(jacobianName).IRight = reshape( ...
        local_normalize(local_i_envelope(rightBasis, n)), 40, 40);
    clusterMaps.(jacobianName).ILeft = reshape( ...
        local_normalize(local_i_envelope(leftBasis, n)), 40, 40);
end

% J_I left eigenvectors have no E component analytically. Remove the
% numerical eigensolver residue instead of normalizing it into a pattern.
clusterMaps.J_I.ELeft = zeros(40, 40);
save(outputFile, 'clusterMaps', 'jacobianNames', 'modeCounts', ...
    'angleDeg', 'contrast', '-v7');
fprintf('Saved Figure 4B maps to %s.\n', outputFile);
fprintf('Cluster mode counts: J_full=%d, J_6=%d, J_I=%d.\n', modeCounts);
end

function [rightBasis, leftBasis] = local_cluster_bases( ...
        matrix, matrixName, modeCount, n, eRows, iRows)
switch matrixName
    case 'J_full'
        active = matrix;
    case 'J_6'
        active = matrix(eRows, eRows);
    case 'J_I'
        active = matrix(iRows, iRows);
    otherwise
        error('Figure4B:Jacobian', 'Unsupported Jacobian %s.', matrixName);
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
    error('Figure4B:Eigs', 'eigs returned flag %d.', flag);
end
values = diag(diagonal);
[~, order] = sortrows([real(values), imag(values)], [-1 -2]);
modes = struct('Vectors', vectors(:, order(1:modeCount)), ...
    'Values', values(order(1:modeCount)));
end

function basis = local_basis(vectors)
vectors = vectors ./ max(vecnorm(vectors), eps);
[basis, triangular] = qr(vectors, 0);
tolerance = max(size(triangular)) * ...
    eps(max(norm(triangular, 2), 1));
rankValue = sum(abs(diag(triangular)) > tolerance);
basis = basis(:, 1:rankValue);
end

function envelope = local_e_envelope(basis, n)
e = 0.6923 * basis(1:n, :) + 0.3077 * basis(n + (1:n), :);
envelope = sqrt(sum(abs(e).^2, 2));
end

function envelope = local_i_envelope(basis, n)
i = basis(2*n + (1:n), :);
envelope = sqrt(sum(abs(i).^2, 2));
end

function values = local_normalize(values)
mapNorm = norm(values(:), 2);
if mapNorm > 1e-12
    values = values / mapNorm;
else
    values = zeros(size(values));
end
end
